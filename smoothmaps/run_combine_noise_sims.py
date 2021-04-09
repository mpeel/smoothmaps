from combine_noise_sims import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits

from smoothmap import *
import os
import matplotlib.pyplot as plt
import scipy.io as io

# Settings
output_resolution = [60.0]#,120.0,240.0]
output_nside = np.asarray([512, 256, 128, 64, 32, 16, 8])
numrealisations = 5
mapnumbers = [0, 1, 2] # II,QQ,UU,QU
rescale = 1.0 # Convert from K to mK - this is applied to the noise map
units_out = 'mKCMB' # Must match rescale!
numres = len(output_resolution)
# directory = '/Users/mpeel/Documents/maps/quijote_202103/reform/'
# directory = '/scratch1/mpeel/maps/'
directory = '/net/nas/proyectos/quijote/validation/Nov2020/noise_simulations/RecommendedSimulations/RecommendedSimulations_allsets_sm1deg/'
# directory = '/share/nas_cbassarc/mpeel/quijote_202103/'
# outdirectory = directory+"../quijote_202103_tqu_noise_v1.0/"
# outdirectory='/Users/mpeel/Desktop/quijote_sims/'
outdirectory='/scratch1/mpeel/quijote_sims/'
os.makedirs(outdirectory, exist_ok=True)

# Window functions
# wf = io.readsav('/Users/mpeel/Documents/maps/quijote_202103/reform/mfi_blconv_wl.sav')
wf = io.readsav('/scratch1/mpeel/mfi_blconv_wl.sav')
# wfq = wf['wl_mfi'].T
wfq = wf['bl'].T
wfq_l = range(0,len(wfq[0][0]))

for i in range(0,numres):
	resolution = "%.2f" % output_resolution[i]

	combine_noise_sims(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI3_11.0_sims_mKCMBunits.fits', prefix='quijote_11GHz_horn3_',postfix='_sm1deg.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[2][0],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,nside_in=512)
	combine_noise_sims(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI3_13.0_sims_mKCMBunits.fits', prefix='quijote_13GHz_horn3_',postfix='_sm1deg.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[2][1],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,nside_in=512)
	combine_noise_sims(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI2_17.0_sims_mKCMBunits.fits', prefix='quijote_17GHz_horn2_',postfix='_sm1deg.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[1][0],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,nside_in=512)
	combine_noise_sims(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI2_19.0_sims_mKCMBunits.fits', prefix='quijote_19GHz_horn2_',postfix='_sm1deg.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[1][1],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,nside_in=512)
	combine_noise_sims(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI4_17.0_sims_mKCMBunits.fits', prefix='quijote_17GHz_horn4_',postfix='_sm1deg.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[3][0],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,nside_in=512)
	combine_noise_sims(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI4_19.0_sims_mKCMBunits.fits', prefix='quijote_19GHz_horn4_',postfix='_sm1deg.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[3][1],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out,nside_in=512)
	# smoothnoisemap(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI3_13.0_2021_mKCMBunits.fits', 'mfi_mar2021_13.0_3.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[2][1],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI2_17.0_2021_mKCMBunits.fits', 'mfi_mar2021_17.0_2.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[1][0],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI2_19.0_2021_mKCMBunits.fits', 'mfi_mar2021_19.0_2.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[1][1],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI4_17.0_2021_mKCMBunits.fits', 'mfi_mar2021_17.0_4.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[3][0],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)
	# smoothnoisemap(directory, outdirectory, resolution+'smoothed_QUIJOTEMFI4_19.0_2021_mKCMBunits.fits', 'mfi_mar2021_19.0_4.fits',mapnumber=mapnumbers,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=wfq[3][1],rescale=rescale,do_intensity=True,do_polarisation=True,units_out=units_out)

# EOF