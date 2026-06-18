import astropy.io.fits as fits
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import os

from smoothnoisemap import *
from smoothmap import *

# Settings
output_resolution = [60.0]
output_nside = np.asarray([512, 256, 128, 64, 32, 16, 8])
numrealisations = 1000
mapnumbers = [3,4, 5, 6] # II,QQ,UU,QU
rescale = 1.0 # Input maps are mK already
units_out = 'mK' # Must match rescale!
numres = len(output_resolution)
directory = '/mnt/DATA1/mpeel/smoothmaps/'
indirectory = directory + "cosmoglobe_dr1/"
outdirectory = directory+"wmap_cosmoglobe_tqu_noise_v1.0_1k/"
os.makedirs(outdirectory, exist_ok=True)

# Read in the beams
beamtf_K = np.loadtxt(indirectory+'wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(indirectory+'wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q1 = np.loadtxt(indirectory+'wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q2 = np.loadtxt(indirectory+'wmap_ampl_bl_Q2_9yr_v5p1.txt',usecols=(1,))
beamtf_V1 = np.loadtxt(indirectory+'wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_V2 = np.loadtxt(indirectory+'wmap_ampl_bl_V2_9yr_v5p1.txt',usecols=(1,))
beamtf_W1 = np.loadtxt(indirectory+'wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))
beamtf_W2 = np.loadtxt(indirectory+'wmap_ampl_bl_W2_9yr_v5p1.txt',usecols=(1,))
beamtf_W3 = np.loadtxt(indirectory+'wmap_ampl_bl_W3_9yr_v5p1.txt',usecols=(1,))
beamtf_W4 = np.loadtxt(indirectory+'wmap_ampl_bl_W4_9yr_v5p1.txt',usecols=(1,))

for i in range(0,numres):
	resolution = "%.2f" % output_resolution[i]

	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmo_22.8_512_2023_mKCMBunits', 'CG_023-WMAP_K_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmo_33.0_512_2023_mKCMBunits', 'CG_030-WMAP_Ka_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmoq1_40.7_512_2023_mKCMBunits', 'CG_040-WMAP_Q1_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmoq2_40.7_512_2023_mKCMBunits', 'CG_040-WMAP_Q2_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmov1_60.7_512_2023_mKCMBunits', 'CG_060-WMAP_V1_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmov2_60.7_512_2023_mKCMBunits', 'CG_060-WMAP_V2_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmoq1_93.5_512_2023_mKCMBunits', 'CG_090-WMAP_W1_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmoq2_93.5_512_2023_mKCMBunits', 'CG_090-WMAP_W2_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmoq3_93.5_512_2023_mKCMBunits', 'CG_090-WMAP_W3_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
	smoothnoisemap(indirectory, outdirectory, resolution+'smoothed_wmap9cosmoq4_93.5_512_2023_mKCMBunits', 'CG_090-WMAP_W4_IQU_n0512_v1.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_K,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,usehealpixfits=True)
