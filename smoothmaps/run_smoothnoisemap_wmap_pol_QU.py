from smoothnoisemap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits

from smoothmap import *
import os
import matplotlib.pyplot as plt

# Settings
output_resolution = [60.0]#,120.0,240.0]
output_nside = [512, 256, 128, 64, 32, 16, 8]
numrealisations = 10
mapnumbers = [0,1,3,2] # II,QQ,UU,QU
hdu=2 # WMAP stores the variance maps in the 2nd header
units_out = 'mK' # MUST MATCH sigma_0 AND sigma_P UNITS BELOW

# directory = '/Users/mpeel/Documents/maps/'
directory = '/scratch1/mpeel/maps/'
outdirectory = directory+"wmap9_tqu_noise_v0.7_TEST/"
os.makedirs(outdirectory, exist_ok=True)
# Read in the beams
beamtf_K = np.loadtxt(directory+'wmap9/wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_V = np.loadtxt(directory+'wmap9/wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_W = np.loadtxt(directory+'wmap9/wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))


for i in range(0,len(output_resolution)):
	resolution = "%.2f" % output_resolution[i]

	sigma_0 = 1.429 # mK
	sigma_P = 1.435 # mK
	smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_22.8_512_2013_mKCMBunits', 'wmap_band_iqumap_r9_9yr_K_v5.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_K,sigma_0=sigma_0,sigma_P=sigma_P,nside=output_nside,hdu=hdu,do_intensity=False,do_polarisation=True,units_out=units_out)

	sigma_0 = 1.466 # mK
	sigma_P = 1.472 # mk
	smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_33.0_512_2013_mKCMBunits', 'wmap_band_iqumap_r9_9yr_Ka_v5.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_Ka,sigma_0=sigma_0,sigma_P=sigma_P,nside=output_nside,hdu=hdu,do_intensity=True,do_polarisation=True,units_out=units_out)

	sigma_0 = 2.188 # mK
	sigma_P = 2.197 # mK
	smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_40.7_512_2013_mKCMBunits', 'wmap_band_iqumap_r9_9yr_Q_v5.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_Q,sigma_0=sigma_0,sigma_P=sigma_P,nside=output_nside,hdu=hdu,do_intensity=True,do_polarisation=True,units_out=units_out)

	sigma_0 = 3.131 # mK
	sigma_P = 3.141 # mK
	smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_60.7_512_2013_mKCMBunits', 'wmap_band_iqumap_r9_9yr_V_v5.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_V,sigma_0=sigma_0,sigma_P=sigma_P,nside=output_nside,hdu=hdu,do_intensity=True,do_polarisation=True,units_out=units_out)

	sigma_0 = 6.544 # mK
	sigma_P = 6.560 # mK
	smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_93.5_512_2013_mKCMBunits', 'wmap_band_iqumap_r9_9yr_W_v5.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_W,sigma_0=sigma_0,sigma_P=sigma_P,nside=output_nside,hdu=hdu,do_intensity=True,do_polarisation=True,units_out=units_out)

