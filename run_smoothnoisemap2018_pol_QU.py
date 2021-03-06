import astropy.io.fits as fits
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import os

from smoothnoisemap import *
from smoothmap import *

def get_hfi_beam(FITSfile):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print(header) # print header
	print(data.names) # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

# Settings
output_resolution = [60.0]#,120.0,240.0]
output_nside = np.asarray([2048, 1024, 512, 256, 128, 64, 32, 16, 8])
numrealisations = 1000
mapnumbers = [4, 7, 9, 8] # II,QQ,UU,QU
rescale = 1000.0 # Convert from K to mK - this is applied to the noise map
units_out = 'mK' # Must match rescale!
numres = len(output_resolution)
# directory = '/Users/mpeel/Documents/maps/'
# directory = '/scratch1/mpeel/maps/'
directory = '/share/nas_cbassarc/mpeel/'
outdirectory = directory+"planck2018_tqu_noise_v1.0/"
os.makedirs(outdirectory, exist_ok=True)

# Read in the beams
beamtf_p30 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',28)
beamtf_p44 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',29)
beamtf_p70 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',30)
beamtf_p100 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_100x100.fits')
beamtf_p143 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_143x143.fits')
beamtf_p217 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_217x217.fits')
beamtf_p353 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_353x353.fits')
beamtf_p545 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_545x545.fits')
beamtf_p857 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_857x857.fits')

for i in range(0,numres):
	resolution = "%.2f" % output_resolution[i]

	# LFI, no bandpass subtraction
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeamnobp_28.4_1024_2018_mKCMBunits', 'LFI_SkyMap_030_1024_R3.00_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p30,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeamnobp_44.1_1024_2018_mKCMBunits', 'LFI_SkyMap_044_1024_R3.00_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p44,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeamnobp_70.4_1024_2018_mKCMBunits', 'LFI_SkyMap_070_1024_R3.00_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p70,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)

	# LFI, with bandpass subtraction
	# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits', 'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p30,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits', 'LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p44,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_70.4_1024_2018_mKCMBunits', 'LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p70,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)

	# # HFI with polarisation
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_100_1024_2018_mKCMBunits', 'HFI_SkyMap_100_2048_R3.01_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p100,rescale=rescale,do_intensity=False,do_polarisation=True,units_out=units_out)
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_143_1024_2018_mKCMBunits', 'HFI_SkyMap_143_2048_R3.01_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p143,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_217_1024_2018_mKCMBunits', 'HFI_SkyMap_217_2048_R3.01_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p217,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_353_1024_2018_mKCMBunits', 'HFI_SkyMap_353-psb_2048_R3.01_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p353,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)

	# Intensity only for 545 and 857GHz
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_545_1024_2018_MJySrunits', 'HFI_SkyMap_545_2048_R3.01_full.fits',mapnumber=[2],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p545,do_intensity=True,do_polarisation=False,units_out='MJySr')
	smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_857_1024_2018_MJySrunits', 'HFI_SkyMap_857_2048_R3.01_full.fits',mapnumber=[2],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p857,do_intensity=True,do_polarisation=False,units_out='MJySr')
