#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Make noise maps from the variance maps, smooth them, and then work out the new variance map
#
# Mike Peel    03 Sep 2017    Start
# Mike Peel    07 Sep 2017    Bug fixes / tidying running order
# Mike Peel    17 Sep 2017    ud_grade to use constant nsides
# Mike Peel    04 Oct 2017    Optimising and debugging. Added multiple nside support.
# Mike Peel    06 Oct 2017    Return to older variance calculation, add rescale param and change output name format
# Mike Peel    05 Jun 2019    v0.4 Add taper, cope with cut sky maps
# Mike Peel    28 Aug 2019    v0.5 Add taper_gauss, normalise, hdu
# Mike Peel    21 Sep 2020    v0.6 Rewrite to use QU covariance matrix
# Mike Peel    28 Sep 2020    v0.7 Correctly convert the WMAP Nobs maps into variance maps

import numpy as np
import numba
import healpy as hp
from smoothmap import smoothmap, conv_nobs_variance_map
import astropy.io.fits as fits
import os
from scipy import optimize

def gaussfit(x, param):
	return hp.gauss_beam(np.radians(param/60.0),300)

# @numba.jit(nopython=True, parallel=True)
def noiserealisation(inputmap, numpixels):
	newmap = np.zeros(numpixels)
	newmap = np.random.normal(scale=1.0, size=numpixels) * inputmap
	return newmap

def precalc_C(Q, U, QU):
	B = np.asarray([[Q,QU],[QU,U]]).T
	print(B)
	# This is test code for finding the covariance array that has a problem
	for arr in B:
		print(arr)
		print(np.linalg.eigvalsh(arr))
		print(np.linalg.cholesky(arr))
	# print(np.min(np.linalg.eigvalsh(B)))
	# print(B[np.argmin(np.linalg.eigvalsh(B))])
	return np.linalg.cholesky(B)

def noiserealisation_QU(C,inputmap_Q,inputmap_U):
	vals = np.random.normal(scale=1.0, size=(len(C),2))
	Q,U = np.einsum('nij,nj->ni', C, vals).T
	return Q, U

def smoothnoisemap(indir, outdir, runname, inputmap, mapnumber=[2], fwhm=0.0, numrealisations=10, sigma_0 = 0.0,sigma_P=0.0, nside=[512], windowfunction = [], rescale=1.0,usehealpixfits=False,taper=False,lmin_taper=350,lmax_taper=600,taper_gauss=False,taper_gauss_sigma=0.0,normalise=True,hdu=1, use_covariance=False, do_intensity=True, do_polarisation=False, units_out='mK'):
	ver = "0.7"

	if (os.path.isfile(indir+"/"+runname+"_actualvariance.fits")):
		print("You already have a file with the output name " + indir+"/"+runname+"_actualvariance.fits" + "! Not going to overwrite it. Move it, or set a new output filename, and try again!")
		# exit()
		return

	# Read in the input map
	print(indir+'/'+inputmap)
	inputfits = fits.open(indir+"/"+inputmap)
	cols = inputfits[hdu].columns
	col_names = cols.names
	nmaps = len(cols)
	maps = []
	if usehealpixfits:
		maps = hp.read_map(indir+inputmap,field=None)
	else:
		for i in range(0,nmaps):
			maps.append(inputfits[hdu].data.field(i))
	nside_in = hp.get_nside(maps)

	# Check to see whether we have nested data, and switch to ring if that is the case.
	if (inputfits[hdu].header['ORDERING'] == 'NESTED'):
		maps = hp.reorder(maps,n2r=True)

	for mapnum in mapnumber:
		maps[mapnum][maps[mapnum]<-1e10] = hp.UNSEEN
		maps[mapnum][~np.isfinite(maps[mapnum])] = hp.UNSEEN
	map_before = maps.copy()

	# If we have Nobs maps, we need to do some preprocessing
	if sigma_0 != 0.0 and do_intensity:
		# Simply convert the intensity map
		maps[mapnumber[0]] = conv_nobs_variance_map(maps[mapnumber[0]], sigma_0)
	if sigma_P != 0.0:
		# Combine the Nobs map into a set of 2x2 matricis
		NobsArr = np.asarray([[maps[mapnumber[1]],maps[mapnumber[3]]],[maps[mapnumber[3]],maps[mapnumber[2]]]]).T
		# Now use the inverse of the matrix, and rescale by sigma_P
		cov = sigma_P * sigma_P * np.linalg.inv(NobsArr)
		# print(cov)
		# Now wrangle it back into the original arrays, now as variances rather than Nobs
		cov = cov.T
		# print(np.shape(cov))
		maps[mapnumber[1]] = cov[0,0]
		maps[mapnumber[2]] = cov[1,1]
		maps[mapnumber[3]] = cov[0,1]

	noisemap = np.zeros((len(mapnumber),len(maps[0])))
	# noisemap = maps.copy()
	i = 0
	for mapnum in mapnumber:
		# if sigma_0 != 0.0:
		# 	# If we have a value for sigma_0, then we have an Nobs map and need to convert it.
		# 	if i > 0 and sigma_P != 0.0:
		# 		# Use the polarisation value
		# 		# Counts in the range -1<0<1 don't make much sense, so zero them (we'll catch the inf later)
		# 		mask = np.zeros(len(maps[mapnum]))
		# 		mask[maps[mapnum]<1.0] = 1.0
		# 		mask[maps[mapnum]<-1.0] = 0.0
		# 		maps[mapnum][mask==1.0] = 0.0
		# 		# Convert
		# 		maps[mapnum] = conv_nobs_variance_map(maps[mapnum], sigma_P)
		# 	else:
		# 		# Use the intensity value
		# 		maps[mapnum] = conv_nobs_variance_map(maps[mapnum], sigma_0)

		# We want to sqrt it to get a noise rms map
		noisemap[i] = np.sqrt(np.abs(maps[mapnum]))
		# Needed for QU, shouldn't make any difference for the others.
		noisemap[i][maps[mapnum]<0] *= -1.0
		noisemap[i] = noisemap[i] * rescale
		# print(np.sum(noisemap[i]<0))
		noisemap[i][map_before[0] == hp.UNSEEN] = 0.0
		noisemap[i][~np.isfinite(noisemap[i])] = 0.0

		# Write the variance map to disk so we can compare to it later.
		cols = []
		outvar = np.square(noisemap[i])
		outvar[noisemap[i]<0] *= -1.0
		cols.append(fits.Column(name='II_cov', format='E', array=outvar))
		cols = fits.ColDefs(cols)
		bin_hdu = fits.BinTableHDU.from_columns(cols)
		bin_hdu.header['ORDERING']='RING'
		bin_hdu.header['POLCONV']='COSMO'
		bin_hdu.header['PIXTYPE']='HEALPIX'
		bin_hdu.header['NSIDE']=nside_in
		bin_hdu.header['TTYPE1'] = 'VARIANCE'
		bin_hdu.header['TUNIT1'] = '('+units_out+")^2"
		bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
		bin_hdu.writeto(outdir+"/"+runname+"_actualvariance_"+str(mapnum)+".fits",overwrite=True)

		# Also save the input nobs map, as a cross-check.
		cols = []
		cols.append(fits.Column(name='II_cov', format='E', array=conv_nobs_variance_map(np.square(noisemap[i]), sigma_0)))
		cols = fits.ColDefs(cols)
		bin_hdu = fits.BinTableHDU.from_columns(cols)
		bin_hdu.header['ORDERING']='RING'
		bin_hdu.header['POLCONV']='COSMO'
		bin_hdu.header['PIXTYPE']='HEALPIX'
		bin_hdu.header['NSIDE']=nside_in
		bin_hdu.header['TTYPE1'] = 'NOBS'
		bin_hdu.header['TUNIT1'] = 'counts'
		bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
		bin_hdu.writeto(outdir+"/"+runname+"_actualnobs_"+str(mapnum)+".fits",overwrite=True)
		i += 1

	# Calculate the window function
	conv_windowfunction = hp.gauss_beam(np.radians(fwhm/60.0),4*nside_in)
	if (windowfunction != []):
		window_len = len(conv_windowfunction)
		beam_len = len(windowfunction)
		if (beam_len > window_len):
			windowfunction  = windowfunction[0:len(conv_windowfunction)]
		else:
			windowfunction = np.pad(windowfunction, (0, window_len - beam_len), 'constant')
		conv_windowfunction[windowfunction!=0] /= windowfunction[windowfunction!=0]
		conv_windowfunction[windowfunction==0] = 0.0
	if normalise:
		conv_windowfunction /= conv_windowfunction[0]
	# If needed, apply a taper
	if taper:
		conv_windowfunction[lmin_taper:lmax_taper] = conv_windowfunction[lmin_taper:lmax_taper] * np.cos((np.pi/2.0)*((np.arange(lmin_taper,lmax_taper)-lmin_taper)/(lmax_taper-lmin_taper)))
		conv_windowfunction[lmax_taper:] = 0.0
	if taper_gauss:
		trip = 0
		val = 0
		# If we haven't been given a FWHM, estimate it from the data at l<300
		if taper_gauss_sigma == 0:
			beam1 = hp.gauss_beam(np.radians(fwhm/60.0),len(conv_windowfunction))
			param_est, cov_x = optimize.curve_fit(gaussfit, range(0,299), windowfunction[0:301], 60.0)
			print(param_est[0])
			taper_gauss_sigma = param_est[0]
		# beam2 = hp.gauss_beam(np.radians(taper_gauss_sigma/60.0),3*nside)
		# plt.plot(windowfunction[0:301])
		# plt.plot(beam2)
		# plt.savefig(outdir+'temp.pdf')
		# exit()
		sigma_final = np.radians(fwhm/60.0)/np.sqrt(8*np.log(2))
		sigma_current = np.radians(taper_gauss_sigma/60.0)/np.sqrt(8*np.log(2))
		for l in range(1,len(conv_windowfunction)):
			if trip == 1:
				conv_windowfunction[l] = val * np.exp(-0.5*(sigma_final**2-sigma_current**2)*l*(l+1))
			elif (conv_windowfunction[l]-conv_windowfunction[l-1]) > 0.0:
				print(l)
				trip = 1
				val = conv_windowfunction[l-1]/np.exp(-0.5*(sigma_final**2-sigma_current**2)*(l-1)*((l-1)+1))
				conv_windowfunction[l] = val * np.exp(-0.5*(sigma_final**2-sigma_current**2)*l*(l+1))


	numpixels = []
	returnmap = []
	num_nside = len(nside)
	for output_nside in nside:
		numpixels.append(hp.nside2npix(output_nside))
		returnmap.append(np.zeros(hp.nside2npix(output_nside)))
	numpixels_orig = len(noisemap[0])


	# Now generate the noise realisations for intensity
	if do_intensity == True:
		noisemap[0][map_before[0] == hp.UNSEEN] = 0.0
		hp.write_map(outdir+"/"+runname+"_noisemap.fits",noisemap[0],overwrite=True)

		for i in range(0,numrealisations):
			if i%10==0:
				print(i)
			# Generate the noise realisation
			newmap = noiserealisation(noisemap[0], numpixels_orig)
			# smooth it
			alms = hp.map2alm(newmap)#,lmax=4*nside_in)
			alms = hp.almxfl(alms, conv_windowfunction)
			newmap = hp.alm2map(alms, nside_in)#,lmax=4*nside_in)
			for j in range(0,num_nside):
				newmap_udgrade = hp.ud_grade(newmap, nside[j], power=0)
				returnmap[j][:] = returnmap[j][:] + np.square(newmap_udgrade)

		for j in range(0,num_nside):
			returnmap[j] = returnmap[j]/(numrealisations-1)
			# returnmap[j][map_before[0] == hp.UNSEEN] = hp.UNSEEN
			
			# All done - now just need to write it to disk.
			cols = []
			cols.append(fits.Column(name='II_cov', format='E', array=returnmap[j]))
			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['NSIDE']=nside[j]
			bin_hdu.header['TTYPE1'] = 'II'
			bin_hdu.header['TUNIT1'] = '('+units_out+")^2"
			bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's smoothnoisemap version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_variance_"+str(nside[j])+".fits",overwrite=True)

			# Also do an Nobs map for a consistency check.
			nobs_map = conv_nobs_variance_map(returnmap[j], sigma_0)
			cols = []
			cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['NSIDE']=nside[j]
			bin_hdu.header['TTYPE1'] = 'Nobs'
			bin_hdu.header['TUNIT1'] = 'count'
			bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_nobs_"+str(nside[j])+".fits",overwrite=True)

		# Do ud_graded versions - no longer used since we generate the different Nsides above
		# num_nside = len(nside)
		# for i in range(0,num_nside):
		# 	# ud_grade it using power=0 (assuming correlated pixels)
		# 	returnmap_ud = hp.ud_grade(returnmap[0], nside[i], power=0)

		# 	# Output the variance map
		# 	cols = []
		# 	cols.append(fits.Column(name='II_cov', format='E', array=returnmap_ud))
		# 	cols = fits.ColDefs(cols)
		# 	bin_hdu = fits.BinTableHDU.from_columns(cols)
		# 	bin_hdu.header['ORDERING']='RING'
		# 	bin_hdu.header['POLCONV']='COSMO'
		# 	bin_hdu.header['PIXTYPE']='HEALPIX'
		# 	bin_hdu.header['NSIDE']=nside[i]
		# 	bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
		# 	bin_hdu.writeto(outdir+"/"+str(nside[i])+"_"+runname+"_variance_udgrade.fits")

		# 	# Also do an Nobs map for a consistency check.
		# 	nobs_map = conv_nobs_variance_map(returnmap_ud, sigma_0)
		# 	cols = []
		# 	cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
		# 	cols = fits.ColDefs(cols)
		# 	bin_hdu = fits.BinTableHDU.from_columns(cols)
		# 	bin_hdu.header['ORDERING']='RING'
		# 	bin_hdu.header['POLCONV']='COSMO'
		# 	bin_hdu.header['PIXTYPE']='HEALPIX'
		# 	bin_hdu.header['NSIDE']=nside[i]
		# 	bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
		# 	bin_hdu.writeto(outdir+"/"+str(nside[i])+"_"+runname+"_nobs_udgrade.fits")

	if do_polarisation == True and len(mapnumber) > 1:
		# ... and now for polarisation
		returnmap_Q = []
		returnmap_U = []
		for output_nside in nside:
			returnmap_Q.append(np.zeros(hp.nside2npix(output_nside)))
			returnmap_U.append(np.zeros(hp.nside2npix(output_nside)))
		hp.write_map(outdir+"/"+runname+"_noisemap_Q.fits",noisemap[1],overwrite=True)
		hp.write_map(outdir+"/"+runname+"_noisemap_U.fits",noisemap[2],overwrite=True)
		hp.write_map(outdir+"/"+runname+"_noisemap_QU.fits",noisemap[3],overwrite=True)

		# Prepare the array
		C = precalc_C(noisemap[1], noisemap[2], noisemap[3])


		for i in range(0,numrealisations):
			if i%10==0:
				print(i)
			# Generate the noise realisation
			newmap_Q, newmap_U = noiserealisation_QU(C,noisemap[1],noisemap[2])
			# smooth it
			alms = hp.map2alm(newmap_Q)#,lmax=4*nside_in)
			alms = hp.almxfl(alms, conv_windowfunction)
			newmap_Q = hp.alm2map(alms, nside_in)#,lmax=4*nside_in)
			alms = hp.map2alm(newmap_U)#,lmax=4*nside_in)
			alms = hp.almxfl(alms, conv_windowfunction)
			newmap_U = hp.alm2map(alms, nside_in)#,lmax=4*nside_in)
			for j in range(0,num_nside):
				newmap_Q_udgrade = hp.ud_grade(newmap_Q, nside[j], power=0)
				newmap_U_udgrade = hp.ud_grade(newmap_U, nside[j], power=0)
				returnmap_Q[j][:] = returnmap_Q[j][:] + np.square(newmap_Q_udgrade)
				returnmap_U[j][:] = returnmap_U[j][:] + np.square(newmap_U_udgrade)

		for j in range(0,num_nside):
			returnmap_Q[j] = returnmap_Q[j]/(numrealisations-1)
			returnmap_U[j] = returnmap_U[j]/(numrealisations-1)
			# returnmap[j][map_before == hp.UNSEEN] = hp.UNSEEN
			
			# All done - now just need to write it to disk.
			cols = []
			cols.append(fits.Column(name='QQ_cov', format='E', array=returnmap_Q[j]))
			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['NSIDE']=nside[j]
			bin_hdu.header['NSIM']=numrealisations
			bin_hdu.header['TTYPE1'] = 'QQ'
			bin_hdu.header['TUNIT1'] = '('+units_out+")^2"
			bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's smoothnoisemap version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_variance_Q_"+str(nside[j])+".fits",overwrite=True)
			cols = []
			cols.append(fits.Column(name='UU_cov', format='E', array=returnmap_U[j]))
			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['NSIDE']=nside[j]
			bin_hdu.header['NSIM']=numrealisations
			bin_hdu.header['TTYPE1'] = 'UU'
			bin_hdu.header['TUNIT1'] = '('+units_out+")^2"
			bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's smoothnoisemap version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_variance_U_"+str(nside[j])+".fits",overwrite=True)

			# Also do an Nobs map for a consistency check.
			nobs_map = conv_nobs_variance_map(returnmap_Q[j], sigma_0)
			cols = []
			cols.append(fits.Column(name='QQ_nobs', format='E', array=nobs_map))
			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['NSIDE']=nside[j]
			bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_nobs_Q_"+str(nside[j])+".fits",overwrite=True)
			nobs_map = conv_nobs_variance_map(returnmap_U[j], sigma_0)
			cols = []
			cols.append(fits.Column(name='UU_nobs', format='E', array=nobs_map))
			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['NSIDE']=nside[j]
			bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_nobs_U_"+str(nside[j])+".fits",overwrite=True)

		# Do ud_graded versions - no longer used since we generate the different Nsides above
		# num_nside = len(nside)
		# for i in range(0,num_nside):
		# 	# ud_grade it using power=0 (assuming correlated pixels)
		# 	returnmap_ud = hp.ud_grade(returnmap[0], nside[i], power=0)

		# 	# Output the variance map
		# 	cols = []
		# 	cols.append(fits.Column(name='II_cov', format='E', array=returnmap_ud))
		# 	cols = fits.ColDefs(cols)
		# 	bin_hdu = fits.BinTableHDU.from_columns(cols)
		# 	bin_hdu.header['ORDERING']='RING'
		# 	bin_hdu.header['POLCONV']='COSMO'
		# 	bin_hdu.header['PIXTYPE']='HEALPIX'
		# 	bin_hdu.header['NSIDE']=nside[i]
		# 	bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
		# 	bin_hdu.writeto(outdir+"/"+str(nside[i])+"_"+runname+"_variance_udgrade.fits")

		# 	# Also do an Nobs map for a consistency check.
		# 	nobs_map = conv_nobs_variance_map(returnmap_ud, sigma_0)
		# 	cols = []
		# 	cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
		# 	cols = fits.ColDefs(cols)
		# 	bin_hdu = fits.BinTableHDU.from_columns(cols)
		# 	bin_hdu.header['ORDERING']='RING'
		# 	bin_hdu.header['POLCONV']='COSMO'
		# 	bin_hdu.header['PIXTYPE']='HEALPIX'
		# 	bin_hdu.header['NSIDE']=nside[i]
		# 	bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
		# 	bin_hdu.writeto(outdir+"/"+str(nside[i])+"_"+runname+"_nobs_udgrade.fits")

	return