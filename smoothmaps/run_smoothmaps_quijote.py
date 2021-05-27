from smoothmap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits
import os
import scipy.io as io
import matplotlib.pyplot as plt

output_resolution = 60.0
output_nside = [512, 256, 128, 64, 32, 16, 8]

# April 2020 recalib3, 4, 5
# directory = '/Users/mpeel/Documents/maps/quijote_202103/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_202103/smooth/'
directory = '/share/nas_cbassarc/mpeel/quijote_202103/'
outdirectory = directory+"../quijote_202103_tqu_v1.5_newwf/"
os.makedirs(outdirectory, exist_ok=True)

wf_new = fits.open(directory+'mfi_btf.fits')
print(wf_new.info())
print(wf_new[1].columns)
print(np.shape(wf_new[1].data[0][0]))
wfq = (wf_new[1].data[0]['BL_CONV'].T)
# exit()

# Window functions
# wf = io.readsav(directory+'mfi_blconv_wl.sav')
# # wfq = wf['wl_mfi'].T
# wfq = wf['bl'].T
# wfq_l = range(0,len(wfq[0][0]))
#
# for i in range(0,4):
# 	plt.plot(wfq[i][0])
# 	plt.plot(wfq_new[i][0])
# 	plt.savefig('test_'+str(i)+'_0.png')
# 	plt.clf()
# 	plt.plot(wfq[i][1])
# 	plt.plot(wfq_new[i][1])
# 	plt.savefig('test_'+str(i)+'_1.png')
# 	plt.clf()
# exit()

prefixes = ['mfi_mar2021','mfi_mar2021_daynight1','mfi_mar2021_daynight2','mfi_mar2021_half1','mfi_mar2021_half2','mfi_mar2021_halfring1','mfi_mar2021_halfring2','mfi_mar2021_period1','mfi_mar2021_period2','mfi_mar2021_period5','mfi_mar2021_period6','mfi_mar2021_pwv1','mfi_mar2021_pwv2']

for k in range(0,len(prefixes)):
	for i in range(0,len(output_nside)):
		ext = prefixes[k].split('_')
		if len(ext) == 3:
			ext = ext[-1]
		else:
			ext = ''
		smoothmap(directory,outdirectory,prefixes[k]+'_11.0_1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI1'+ext+'_11.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[0][0],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)
		smoothmap(directory,outdirectory,prefixes[k]+'_13.0_1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI1'+ext+'_13.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[0][1],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)
		smoothmap(directory,outdirectory,prefixes[k]+'_11.0_3.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI3'+ext+'_11.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[2][0],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)
		smoothmap(directory,outdirectory,prefixes[k]+'_13.0_3.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI3'+ext+'_13.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[2][1],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)
		smoothmap(directory,outdirectory,prefixes[k]+'_17.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI2'+ext+'_17.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[1][0],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)
		smoothmap(directory,outdirectory,prefixes[k]+'_19.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI2'+ext+'_19.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[1][1],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)
		smoothmap(directory,outdirectory,prefixes[k]+'_17.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI4'+ext+'_17.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[3][0],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)
		smoothmap(directory,outdirectory,prefixes[k]+'_19.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_QUIJOTEMFI4'+ext+'_19.0_2021_mKCMBunits.fits', output_resolution,windowfunction=wfq[3][1],nside_out=output_nside[i],units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True)

# # EOF
