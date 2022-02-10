#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Combine noise simulations to get the smoothed variance maps
#
# Mike Peel    09 Apr 2021    v1.0 Fork from smoothnoisemap
# NB: Known bug, the rescale factor does not seem to work!

import numpy as np
import numba
import healpy as hp
from smoothmap import smoothmap, conv_nobs_variance_map
import astropy.io.fits as fits
import os
from scipy import optimize

def combine_noise_sims(indir, outdir, runname, prefix='',postfix='', mapnumber=[], fwhm=0.0, numrealisations=10, sigma_0 = 0.0,sigma_P=0.0, nside=[512], windowfunction = [], rescale=1.0,usehealpixfits=False,taper=False,lmin_taper=350,lmax_taper=600,taper_gauss=False,taper_gauss_sigma=0.0,normalise=True,hdu=1, use_covariance=False, do_intensity=True, do_polarisation=False, units_out='mK',do_smoothing=True,use_precomputed_wf=False,nside_in=512):
	ver = "1.0"

	# Calculate the window function
	if use_precomputed_wf:
		print('Using precomputed window function')
		conv_windowfunction = windowfunction
		conv_windowfunction = np.pad(conv_windowfunction, (0, 3*nside - len(conv_windowfunction)), 'constant')
	else:
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

	returnmap_Q = []
	returnmap_U = []
	returnmap_QU = []
	for output_nside in nside:
		returnmap_Q.append(np.zeros(hp.nside2npix(output_nside)))
		returnmap_U.append(np.zeros(hp.nside2npix(output_nside)))
		returnmap_QU.append(np.zeros(hp.nside2npix(output_nside)))


	for i in range(0,numrealisations):
		if i%10==0:
			print(i)

		# Read in the input map
		print(indir+'/'+prefix+str(i+1).zfill(4)+postfix)
		maps = []
		if usehealpixfits:
			maps = hp.read_map(indir+prefix+str(i+1).zfill(4)+postfix,field=None)
		else:
			inputfits = fits.open(indir+"/"+prefix+str(i+1).zfill(4)+postfix)
			cols = inputfits[hdu].columns
			col_names = cols.names
			nmaps = len(cols)
			for i in range(0,nmaps):
				maps.append(inputfits[hdu].data.field(i))
			# Check to see whether we have nested data, and switch to ring if that is the case.
			try:
				if (inputfits[hdu].header['ORDERING'] == 'NESTED'):
					maps = hp.reorder(maps,n2r=True)
			except:
				null = 0

		alms = hp.map2alm([maps[0], maps[1], maps[2]], pol=True)
		alms = hp.sphtfunc.smoothalm(alms, beam_window=conv_windowfunction,pol=True)
		newmaps = hp.alm2map(alms,nside_in,pol=True)

		for j in range(0,num_nside):
			if do_intensity:
				newmap_udgrade = hp.ud_grade(newmaps[0], nside[j], power=0)
				returnmap[j][:] = returnmap[j][:] + np.square(newmap_udgrade)
			if do_polarisation:
				newmap_Q_udgrade = hp.ud_grade(newmaps[1], nside[j], power=0)
				newmap_U_udgrade = hp.ud_grade(newmaps[2], nside[j], power=0)
				returnmap_Q[j][:] = returnmap_Q[j][:] + np.square(newmap_Q_udgrade)
				returnmap_U[j][:] = returnmap_U[j][:] + np.square(newmap_U_udgrade)
				returnmap_QU[j][:] = returnmap_QU[j][:] + newmap_Q_udgrade*newmap_U_udgrade

	for j in range(0,num_nside):
		if do_intensity:
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
			bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's combine_noise_sims version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_variance_"+str(nside[j])+".fits",overwrite=True)

		if do_polarisation:
			returnmap_Q[j] = returnmap_Q[j]/(numrealisations-1)
			returnmap_U[j] = returnmap_U[j]/(numrealisations-1)
			returnmap_QU[j] = returnmap_QU[j]/(numrealisations-1)
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
			bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's combine_noise_sims version "+ver +"."
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
			bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's combine_noise_sims version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_variance_U_"+str(nside[j])+".fits",overwrite=True)
			cols = []
			cols.append(fits.Column(name='QU_cov', format='E', array=returnmap_QU[j]))
			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['NSIDE']=nside[j]
			bin_hdu.header['NSIM']=numrealisations
			bin_hdu.header['TTYPE1'] = 'QU'
			bin_hdu.header['TUNIT1'] = '('+units_out+")^2"
			bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's combine_noise_sims version "+ver +"."
			bin_hdu.writeto(outdir+"/"+runname+"_variance_QU_"+str(nside[j])+".fits",overwrite=True)

	return
