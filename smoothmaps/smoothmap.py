#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Read in a map, smooth it, and write it out
#
# History:
# v0.1  Mike Peel   10-Jan-2016   Initial version.
# v0.2  Mike Peel   8-Mar-2016    Start to generalise to cover other numbers of maps
# v0.3  Mike Peel   6-Sep-2016    Switch to using pyfits; use alms and window functions to smooth; add unit conversion functionality
# v0.4  Mike Peel   23-Sep-2016   Carry fits headers through to the output file. Make fwhm_arcmin optional so that the function can be used solely to ud_grade maps.
# v0.5  Mike Peel   23-Sep-2016   Adding calc_variance_windowfunction (port of Paddy Leahy's IDL code) to properly smooth variance maps.
# v0.6  Mike Peel   26-Sep-2016   Support Nobs maps, and conversion to variance maps, plus debugging/tidying.
# v0.7  Mike Peel   27-Jan-2017   Bug fix - ud_grading the variance maps should include a factor of (nside_orig/nside_new)^2
# V0.8  Adam Barr   Various       Various tweaks
# v0.9  Mike Peel   17 Sep 2017   Tweaks, add rescale parameter
# v0.9a Mike Peel   21 Sep 2017   More bug fixes
# v0.9b Mike Peel   22 Sep 2017   Add nosmooth and outputmaps parameters
# v0.9c Mike Peel   29 Sep 2017   More bug fixes
# v0.9d Mike Peel   01 Oct 2017   Input/output dir improvements. Fix bug in beam window functions.
# v1.0  Mike Peel   06 Oct 2017   Version 1.0. Add option to append maps (e.g. MC output)
# v1.1  Mike Peel   09 Oct 2017   Add CMB subtraction
# v1.2  Mike Peel   04 Jul 2019   Multiple options for window function tapering, min/max map values
# v1.3  Mike Peel   28 Aug 2019   Gaussian taper, options for normalising wf and unseen vs. 0 in map
# v1.4  Mike Peel   24 Jul 2020   Add an option to not smooth the variance maps (but save them in the output anyway)
# v1.4a Mike Peel   27 Jul 2020   Tweak to only check for NESTED when not using usehealpixfits
# v1.4b Mike Peel   07 Oct 2020   Tweak to use 'usehealpixfits' for 'subtractmap'
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import scipy as sp
import math as m

# Work-around for no X display
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from astrocode.fitspectrum.spectra import *
import astropy.io.fits as fits
from scipy import special
import os.path
from scipy import optimize

def gaussfit(x, param):
	return hp.gauss_beam(np.radians(param/60.0),300)


def smoothmap(indir, outdir, inputfile, outputfile, fwhm_arcmin=-1, nside_out=0,maxnummaps=-1, frequency=100.0, units_in='',units_out='', windowfunction = [],nobs_out=False,variance_out=True, sigma_0 = -1, sigma_0_unit='', rescale=1.0, nosmooth=[], outputmaps=[],appendmap='',appendmapname='',appendmapunit='',subtractmap='',subtractmap_units='',usehealpixfits=False,taper=False,lmin_taper=350,lmax_taper=600, cap_one=False, cap_oneall=False,minmapvalue=0,maxmapvalue=0,minmaxmaps=[0],taper_gauss=False,taper_gauss_sigma=0.0,normalise=True,useunseen=False,smoothvariance=False):
	ver = "1.4"

	if (os.path.isfile(outdir+outputfile)):
		print("You already have a file with the output name " + outdir+outputfile + "! Not going to overwrite it. Move it, or set a new output filename, and try again!")
		# exit()
		return

	# Check to see if we have a sigma_0 value to use when converting from Nobs maps and back.
	no_sigma_0 = False
	if (sigma_0 == -1):
		no_sigma_0 = True

	# Read in the fits map, and put it into the format Healpix expects
	# if usehealpixfits:
	# 	maps, newheader = hp.read_map(indir+inputfile,field=None,h=True)
	# 	nmaps = len(maps)
	# 	print newheader
	# else:
	inputfits = fits.open(indir+inputfile)
	print(inputfits.info())
	cols = inputfits[1].columns
	col_names = cols.names
	nmaps = len(cols)
	maps = []
	if usehealpixfits:
		maps = hp.read_map(indir+inputfile,field=None)
	else:
		for i in range(0,nmaps):
			maps.append(inputfits[1].data.field(i))
			print(len(maps[i]))
		print(len(maps[0]))
		# Check to see whether we have nested data, and switch to ring if that is the case.
		if (inputfits[1].header['ORDERING'] == 'NESTED'):
			maps = hp.reorder(maps,n2r=True)
	newheader = inputfits[1].header.copy(strip=False)
	inputfits.close()

	# Do some cleanup of the new header to avoid issues later
	try:
		print(newheader['TUNIT1'])
	except:
		for i in range(0,nmaps):
			newheader['TUNIT'+str(i+1)] = 'None'


	# Crop to just have the maps we want to output
	nmaps_orig = nmaps
	if maxnummaps != -1:
		nmaps = maxnummaps
	if outputmaps == []:
		outputmaps = range(0,nmaps)
	maps = [maps[i] for i in outputmaps]
	noutputmaps = len(outputmaps)
	# print newheader
	for i in range(0,nmaps_orig):
		if i < noutputmaps:
			newheader['TTYPE'+str(i+1)] = newheader['TTYPE'+str(outputmaps[i]+1)]
			newheader['TFORM'+str(i+1)] = newheader['TFORM'+str(outputmaps[i]+1)]
			newheader['TUNIT'+str(i+1)] = newheader['TUNIT'+str(outputmaps[i]+1)]
		else:
			del newheader['TTYPE'+str(i+1)]
			del newheader['TFORM'+str(i+1)]
			del newheader['TUNIT'+str(i+1)]
	newheader['TFIELDS'] = noutputmaps
	newheader['NAXIS1'] = noutputmaps*4
	nmaps = noutputmaps
	print(newheader)

	# Calculate the unit conversion factor
	const = get_spectrum_constants()

	pix_area = hp.nside2pixarea(nside_out)

	nside = hp.get_nside(maps)
	if (fwhm_arcmin != -1):
		conv_windowfunction = hp.gauss_beam(np.radians(fwhm_arcmin/60.0),3*nside)
		if (windowfunction != []):
			window_len = len(conv_windowfunction)
			beam_len = len(windowfunction)

			if (beam_len > window_len):
				windowfunction  = windowfunction[0:len(conv_windowfunction)]
			else:
				# We want to pad the window function rather than crop the convolution
				# conv_windowfunction = conv_windowfunction[0:len(windowfunction)]
				windowfunction = np.pad(windowfunction, (0, window_len - beam_len), 'constant')

			conv_windowfunction[windowfunction!=0] /= windowfunction[windowfunction!=0]
			conv_windowfunction[windowfunction==0] = 0.0
		# Normalise window function
		if normalise:
			conv_windowfunction /= conv_windowfunction[0]
		
		conv_windowfunction_before = conv_windowfunction.copy()
		# If needed, apply a taper
		if taper:
			conv_windowfunction[lmin_taper:lmax_taper] = conv_windowfunction[lmin_taper:lmax_taper] * np.cos((np.pi/2.0)*((np.arange(lmin_taper,lmax_taper)-lmin_taper)/(lmax_taper-lmin_taper)))
			conv_windowfunction[lmax_taper:] = 0.0
		if cap_one:
			conv_windowfunction[conv_windowfunction >= 1.0] = 1.0
		if cap_oneall:
			trip = 0
			for l in range(0,len(conv_windowfunction)):
				if trip == 1:
					conv_windowfunction[l] = 1.0
				elif conv_windowfunction[l] >= 1.0:
					trip = 1
		if taper_gauss:
			trip = 0
			val = 0
			# If we haven't been given a FWHM, estimate it from the data at l<300
			if taper_gauss_sigma == 0:
				beam1 = hp.gauss_beam(np.radians(fwhm_arcmin/60.0),3*nside)
				param_est, cov_x = optimize.curve_fit(gaussfit, range(0,299), windowfunction[0:301], 60.0)
				print(param_est[0])
				taper_gauss_sigma = param_est[0]
			# beam2 = hp.gauss_beam(np.radians(taper_gauss_sigma/60.0),3*nside)
			# plt.plot(windowfunction[0:301])
			# plt.plot(beam2)
			# plt.savefig(outdir+'temp.pdf')
			# exit()
			sigma_final = np.radians(fwhm_arcmin/60.0)/np.sqrt(8*np.log(2))
			sigma_current = np.radians(taper_gauss_sigma/60.0)/np.sqrt(8*np.log(2))
			for l in range(1,len(conv_windowfunction)):
				if trip == 1:
					conv_windowfunction[l] = val * np.exp(-0.5*(sigma_final**2-sigma_current**2)*l*(l+1))
					# conv_windowfunction[l] = val * (beam1[l]/beam2[l])
				elif (conv_windowfunction[l]-conv_windowfunction[l-1]) > 0.0:
					print(l)
					trip = 1
					val = conv_windowfunction[l-1]/np.exp(-0.5*(sigma_final**2-sigma_current**2)*(l-1)*((l-1)+1))
					conv_windowfunction[l] = val * np.exp(-0.5*(sigma_final**2-sigma_current**2)*l*(l+1))
					# val = conv_windowfunction[l-1]/(beam1[l-1]/beam2[l-1])
					# conv_windowfunction[l] = val * (beam1[l]/beam2[l])


		# plt.xscale('linear')
		# plt.yscale('log')
		# plt.plot(windowfunction,label='Window function')
		# plt.plot(hp.gauss_beam(np.radians(fwhm_arcmin/60.0),3*nside),label='Gaussian')
		# plt.legend()
		# plt.savefig(outdir+'wf_orig_'+outputfile+'.pdf')
		# plt.clf()

		# plt.xscale('linear')
		# plt.yscale('linear')
		# plt.ylim(0.0,1.0)
		# plt.plot(conv_windowfunction,label='Window function')
		# plt.plot(conv_windowfunction_before,label='Window function before')
		# plt.plot(hp.gauss_beam(np.radians(fwhm_arcmin/60.0),len(windowfunction)-1)/windowfunction,label='Comparison')
		# plt.legend()
		# plt.savefig(outdir+'wf_lin_'+outputfile+'.pdf')
		# plt.yscale('log')
		# plt.ylim(1e-8,1e4)
		# plt.savefig(outdir+'wf_log_'+outputfile+'.pdf')
		# plt.clf()
		# return

		# Check whether we'll need to smooth variances too.
		if smoothvariance != False:
			test = False
			for i in range(0,nmaps):
				if ('cov' in newheader['TTYPE'+str(i+1)]) or ('N_OBS' in newheader['TTYPE'+str(i+1)]):
					test = True
			if test:
				print('Covariance maps detected. Calculating variance window function (this may take a short while)')
				conv_windowfunction_variance = conv_windowfunction.copy()
				# conv_windowfunction_variance = calc_variance_windowfunction(conv_windowfunction)
				# conv_windowfunction_variance /= conv_windowfunction_variance[0]
				print(conv_windowfunction_variance[0])
				print('Done! Onwards...')

				plt.xscale('log')
				plt.yscale('log')
				plt.plot(conv_windowfunction,label='Window function')
				plt.plot(conv_windowfunction_variance,label='Variance window function')
				plt.legend()
				plt.savefig(outdir+'wf_'+outputfile+'.png')
				plt.xscale('log')
				plt.yscale('linear')
				plt.plot(conv_windowfunction,label='Window function')
				plt.plot(conv_windowfunction_variance,label='Variance window function')
				plt.legend()
				plt.savefig(outdir+'wf_lin_'+outputfile+'.png')

	# Do the smoothing
	print("Smoothing the maps")
	smoothed_map = maps
	for i in range(0,nmaps):
		print('map ' + str(i))
		print(len(maps[i]))
		map_before = maps[i][:].copy()
		if useunseen == True:
			maps[i][~np.isfinite(maps[i][:])] = hp.UNSEEN
			# See if we want to cut based on min/max map values
			if minmapvalue != maxmapvalue:
				if i in minmaxmaps:
					maps[i][maps[i][:] < minmapvalue] = hp.UNSEEN
					maps[i][maps[i][:] > maxmapvalue] = hp.UNSEEN
		else:
			maps[i][maps[i][:] == hp.UNSEEN] = 0.0
			maps[i][~np.isfinite(maps[i][:])] = 0.0
			# See if we want to cut based on min/max map values
			if minmapvalue != maxmapvalue:
				if i in minmaxmaps:
					maps[i][maps[i][:] < minmapvalue] = 0.0
					maps[i][maps[i][:] > maxmapvalue] = 0.0

		# Check that we actually want to do smoothing, as opposed to udgrading. Also check to see if this is in the list of maps to not smooth
		if fwhm_arcmin != -1 and (i not in nosmooth):
			if 'N_OBS' in newheader['TTYPE'+str(i+1)]:
				print('Column '+str(i)+' is an N_OBS map ('+newheader['TUNIT'+str(i+1)]+') - converting to variance map.')
				print(np.sum(maps[i]))
				print(np.median(maps[i]))
				maps[i] = conv_nobs_variance_map(maps[i], sigma_0)
				print(np.sum(maps[i]))
				print(np.median(maps[i]))
				if (nobs_out == False and no_sigma_0 == False):
					# We don't want to convert back later.
					print('test')
					newheader['TTYPE'+str(i+1)] = 'II_cov'
					newheader['TUNIT'+str(i+1)] = '('+sigma_0_unit+')^2'

			# Calculate the alm's, multiply them by the window function, and convert back to the map
			if ('cov' in newheader['TTYPE'+str(i+1)]) or ('N_OBS' in newheader['TTYPE'+str(i+1)]):
				if smoothvariance != False:
					print('Column '+str(i)+' is a covariance matrix ('+newheader['TUNIT'+str(i+1)]+') - smoothing appropriately.')
					alms = hp.map2alm(maps[i])
					alms = hp.almxfl(alms, conv_windowfunction_variance)
					newmap = hp.alm2map(alms, nside,verbose=False)
				else:
					newmap = maps[i].copy()
			else:
				alms = hp.map2alm(maps[i])
				alms = hp.almxfl(alms, conv_windowfunction)
				newmap = hp.alm2map(alms, nside,verbose=False)
			smoothed_map[i] = newmap
			print(np.sum(smoothed_map[i]))
			print(np.median(smoothed_map[i]))
			smoothed_map[i][map_before[:] == hp.UNSEEN] = hp.UNSEEN

			if ('N_OBS' in newheader['TTYPE'+str(i+1)]) and (nobs_out or no_sigma_0):
				print('You\'ve either asked for an N_OBS map to be returned, or not set sigma_0, so you will get an N_OBS map returned in your data!')
				print(np.sum(smoothed_map[i]))
				smoothed_map[i] = conv_nobs_variance_map(smoothed_map[i], sigma_0)
				print(np.sum(smoothed_map[i]))
				newheader['TTYPE'+str(i+1)] = 'N_OBS'
	maps = 0
	newmap = 0

	# Do the ud_grading
	print("ud_grading the maps (if needed)")
	nobs_sum = 0
	if (nside_out == 0):
		nside_out = nside
	else:
		for i in range (0,nmaps):
			if 'N_OBS' in newheader['TTYPE'+str(i+1)]:
				nobs_sum = np.sum(smoothed_map[i])

			# Check to see which type of map we have, and adjust the factor of (nside/nside_out)^power appropriately
			power = 0.0
			if ('cov' in newheader['TTYPE'+str(i+1)]):
				power = 2.0
			elif ('N_OBS' in newheader['TTYPE'+str(i+1)]) or ('Hits' in newheader['TTYPE'+str(i+1)]):
				power = -2.0
			print(power)
			smoothed_map[i] = hp.ud_grade(smoothed_map[i], nside_out=nside_out, power=power)

			if 'N_OBS' in newheader['TTYPE'+str(i+1)]:
				null = 0
			else:
				# If we don't have an N_OBS map, then we might want to convert the units.
				if (units_out != ''):
					if (units_in == ''):
						unit = newheader['TUNIT'+str(i+1)]
						unit = unit.strip()
					else:
						# Assume the user is right to have specified different input units from what is in the file.
						unit = units_in

					newheader['TUNIT'+str(i+1)] = units_out
					power = 1.0
					if ('^2' in unit):
						power = 2.0
						unit = unit.replace(")^2",'').replace('(','')
						unit = unit.replace("^2",'')
						newheader['TUNIT'+str(i+1)] = '('+units_out+")^2"
					print(unit + " " + str(power))
					conversion = convertunits(const, unit, units_out, frequency, pix_area)
					print(conversion)
					smoothed_map[i] = smoothed_map[i] * conversion**power

	# All done - now just need to write it to disk.
	print("Writing maps to disk: " + outdir+outputfile)
	cols = []
	for i in range(0,nmaps):
		if ('cov' in newheader['TTYPE'+str(i+1)]):
			smoothed_map[i] = smoothed_map[i] * rescale**2
		else:
			smoothed_map[i] = smoothed_map[i] * rescale

	# If we want to subtract another map (e.g., a CMB map) then we need to read it in, check nside and units, and then subtract.
	if subtractmap != '':
		print('Subtracting CMB map ' + subtractmap)
		if usehealpixfits:
			sub_maps = hp.read_map(indir+subtractmap)
		else:
			sub_inputfits = fits.open(indir+subtractmap)
			sub_cols = sub_inputfits[1].columns
			sub_col_names = sub_cols.names
			sub_nmaps = len(sub_cols)
			sub_nmaps_orig = sub_nmaps
			sub_maps = []
			for i in range(0,sub_nmaps):
				sub_maps.append(sub_inputfits[1].data.field(i))
			# Check to see whether we have nested data, and switch to ring if that is the case.
			if (sub_inputfits[1].header['ORDERING'] == 'NESTED'):
				sub_maps = hp.reorder(sub_maps,n2r=True)
			# We want a maximum of one map to subtract (could be extended in the future)
			sub_maps = sub_maps[0]
			sub_inputfits.close()
		# We want to use the same output Nside - we'll assume they have the same resolution.
		sub_maps = hp.ud_grade(sub_maps, nside_out=nside_out)
		# Calculate a rescaling factor if needed
		conversion = 1.0
		if subtractmap_units != '':
			conversion = convertunits(const, subtractmap_units, units_out, frequency, pix_area)
		# ... and do the subtraction
		smoothed_map[0] = smoothed_map[0] - sub_maps * conversion
	
	for i in range(0,nmaps):
		cols.append(fits.Column(name=col_names[i], format='E', array=smoothed_map[i]))

	if appendmap != '':
		addmap = hp.read_map(appendmap)
		cols.append(fits.Column(name=appendmapname, format='E', array=addmap))
		newheader['TTYPE'+str(nmaps+1)] = appendmapname
		newheader['TFORM'+str(nmaps+1)] = 'E'
		newheader['TTYPE'+str(nmaps+1)] = appendmapname
		newheader['TFIELDS'] = newheader['TFIELDS']+1
		newheader['NAXIS1'] = newheader['TFIELDS']*4


	cols = fits.ColDefs(cols)
	bin_hdu = fits.BinTableHDU.from_columns(cols)
	# bin_hdu = fits.new_table(cols)
	# print newheader
	bin_hdu.header = newheader
	# print bin_hdu.header
	bin_hdu.header['ORDERING']='RING'
	bin_hdu.header['POLCONV']='COSMO'
	bin_hdu.header['PIXTYPE']='HEALPIX'
	bin_hdu.header['NSIDE']=nside_out
	bin_hdu.header['COMMENT']="Smoothed using Mike Peel's smoothmap.py version "+ver +" modified by Adam Barr"
	print(bin_hdu.header)
	
	bin_hdu.writeto(outdir+outputfile)

	return

def conv_nobs_variance_map(inputmap, sigma_0):
	newmap = sigma_0**2 / inputmap
	return newmap

# Code to replicate IDL's INT_TABULATED function
# From http://stackoverflow.com/questions/14345001/idls-int-tabulate-scipy-equivalent
def int_tabulated(x, f, p=5) :
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = sp.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in xrange(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret

def calc_variance_windowfunction(conv_windowfunction):
	# Calculate the window function for variance maps.
	# Based on IDL code by Paddy Leahy, 'write_cvbl_planck3.pro'

	const = get_spectrum_constants()
	nbl = len(conv_windowfunction)
	ll = np.arange(0,nbl,1)
	
	# Choose scale to match size of beam. First find half-power point
	lhalf = [ n for n,i in enumerate(conv_windowfunction) if i<0.5 ][0]
	
	# Calculate beam out to about 40 * half-power radius, roughly (note that
	# this gives 100 points out to the half-power point).
	numelements = 4000
	rad = np.arange(0,numelements,1)*10.0*const['pi']/(float(lhalf)*float(numelements))
	
	x = np.cos(rad)
	sinrad = np.sin(rad)

	lgndr = np.zeros((numelements,nbl))
	for i in range(0,nbl):
		lgndr[:,i] = special.lpmv(0, i, x)
	
	# Generate radial profile of convolving beam:
	conva = np.zeros(numelements)
	for j in range(0,numelements):
		conva[j] = np.sum((ll+0.5)*conv_windowfunction*lgndr[j,:])

	conva = conva / (2.0*const['pi'])
	print('Peak of convolving beam is ' + str(conva[0]) + " (check: " + str(np.max(conva)) + ")")

	# Square convolving beam and convert back to window function
	mult = sinrad*conva**2
	cvbl = np.zeros(nbl)
	print(nbl)
	for l in range(0,nbl):
		cvbl[l] = int_tabulated(rad,mult*lgndr[:,l])

	# Put in 2pi normalization factor:
	cvbl = 2.0*const['pi']*cvbl

	print('Max in the window function is ' + str(cvbl[0]) + " (check: " + str(np.max(cvbl)) + ")")

	return cvbl


# Read in a Planck (HFI or LFI) beam
def get_beam(FITSfile,hdu=0):
	# fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, hdu, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	# print(header) # print header
	# print(data.names) # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

# EOF
