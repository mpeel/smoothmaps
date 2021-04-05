from smoothnoisemap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits

from smoothmap import *
import os
import matplotlib.pyplot as plt

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
mapnumbers_bpcorrect = [3,4,5,6] # Because the bp-corrected ones are different
rescale = 1000.0 # Convert from K to mK - this is applied to the noise map
units_out = 'mK' # Must match rescale!
numres = len(output_resolution)
# directory = '/Users/mpeel/Documents/maps/'
# directory = '/scratch1/mpeel/maps/'
directory = '/share/nas_cbassarc/mpeel/'
outdirectory = directory+"planck2015_tqu_noise_v1.0/"
os.makedirs(outdirectory, exist_ok=True)

# Read in the beams
beamtf_p30 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',28)
beamtf_p44 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',29)
beamtf_p70 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',30)

HFIbeams = fits.open(directory+'planck2015/HFI_RIMO_Beams-100pc_R2.00.fits')
beamtf_p100 = HFIbeams[3].data[0][0]
beamtf_p143 = HFIbeams[4].data[0][0]
beamtf_p217 = HFIbeams[5].data[0][0]
beamtf_p353 = HFIbeams[6].data[0][0]
beamtf_p545 = HFIbeams[7].data[0][0]
beamtf_p857 = HFIbeams[8].data[0][0]

for i in range(0,numres):
	resolution = "%.2f" % output_resolution[i]

	# LFI, no bandpass subtraction
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_28.4_1024_2015_mKCMBunits', 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p30,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_44.1_1024_2015_mKCMBunits', 'LFI_SkyMap_044_1024_R2.01_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=1024],windowfunction=beamtf_p44,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_70.4_1024_2015_mKCMBunits', 'LFI_SkyMap_070_2048_R2.01_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p70,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)

	# LFI, with bandpass subtraction
	# Note that the signal maps are smoothed - assuming the variance maps aren't!
	# ... can't assume that, the variance maps are also at 1°!
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeambpcorr_28.4_1024_2015_mKCMBunits', 'LFI_SkyMap_030-BPassCorrected_0256_R2.01_full.fits',mapnumber=mapnumbers_bpcorrect,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=256],windowfunction=beamtf_p30,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeambpcorr_44.1_1024_2015_mKCMBunits', 'LFI_SkyMap_044-BPassCorrected_0256_R2.01_full.fits',mapnumber=mapnumbers_bpcorrect,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=256],windowfunction=beamtf_p44,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeambpcorr_70.4_1024_2015_mKCMBunits', 'LFI_SkyMap_070-BPassCorrected_0256_R2.01_full.fits',mapnumber=mapnumbers_bpcorrect,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=256],windowfunction=beamtf_p70,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)

	# # HFI with polarisation
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_100_1024_2015_mKCMBunits', 'HFI_SkyMap_100_2048_R2.02_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p100,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_143_1024_2015_mKCMBunits', 'HFI_SkyMap_143_2048_R2.02_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p143,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_217_1024_2015_mKCMBunits', 'HFI_SkyMap_217_2048_R2.02_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p217,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_353_1024_2015_mKCMBunits', 'HFI_SkyMap_353_2048_R2.02_full.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p353,rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)

	# # Intensity only for 545 and 857GHz
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_545_1024_2015_MJySrunits', 'HFI_SkyMap_545_2048_R2.02_full.fits',mapnumber=[2],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p545,do_intensity=True,do_polarisation=False,units_out='MJySr')
	# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_857_1024_2015_MJySrunits', 'HFI_SkyMap_857_2048_R2.02_full.fits',mapnumber=[2],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p857,do_intensity=True,do_polarisation=False,units_out='MJySr')

	# The bandpass correction variances - already at 1°, but we want to get ud_graded versions
	smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2bpasscorrection_30_mKCMBunits', 'LFI_CorrMap_030-BPassCorr_0256_R2.01_full.fits',mapnumber=mapnumbers_bpcorrect,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=256],windowfunction=beamtf_p30,do_intensity=False,do_polarisation=True,units_out=units_out,do_smoothing=False,usehealpixfits=True)
	smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2bpasscorrection_44_mKCMBunits', 'LFI_CorrMap_044-BPassCorr_0256_R2.01_full.fits',mapnumber=mapnumbers_bpcorrect,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=256],windowfunction=beamtf_p44,do_intensity=False,do_polarisation=True,units_out=units_out,do_smoothing=False,usehealpixfits=True)
	smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2bpasscorrection_70_mKCMBunits', 'LFI_CorrMap_070-BPassCorr_0256_R2.01_full.fits',mapnumber=mapnumbers_bpcorrect,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside[output_nside<=256],windowfunction=beamtf_p70,do_intensity=False,do_polarisation=True,units_out=units_out,do_smoothing=False,usehealpixfits=True)
