import astropy.io.fits as fits
import healpy as hp
import numpy as np
import os

from smoothmap import *

output_resolution = 60.0
output_nside = [512, 256, 128, 64, 32, 16, 8]
smoothvariance = False

directory = '/mnt/DATA1/mpeel/smoothmaps/'
indirectory = directory + "cosmoglobe_dr1/"
outdirectory = directory + "wmap9_cosmoglobe_tqu_v1.5/"
os.makedirs(outdirectory, exist_ok=True)

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

numnside = len(output_nside)
for i in range(0,numnside):
	subtractmaps = ['']
	subtractmaps_name = ['']
	numsubtract = len(subtractmaps)
	for j in range(0,numsubtract):

		if output_nside[i] <= 512:
			# NB: these only use sigma_0 and Nobs for the intensity maps!
			smoothmap(indirectory,outdirectory,'CG_023-WMAP_K_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmo'+subtractmaps_name[j]+'_22.8_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_K,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_030-WMAP_Ka_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmo'+subtractmaps_name[j]+'_33.0_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_Ka,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_040-WMAP_Q1_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmoq1'+subtractmaps_name[j]+'_40.7_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_Q1,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_040-WMAP_Q2_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmoq2'+subtractmaps_name[j]+'_40.7_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_Q2,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_060-WMAP_V1_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmov1'+subtractmaps_name[j]+'_60.7_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_V1,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_060-WMAP_V2_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmov2'+subtractmaps_name[j]+'_60.7_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_V2,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_090-WMAP_W1_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmow1'+subtractmaps_name[j]+'_93.5_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_W1,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_090-WMAP_W2_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmow2'+subtractmaps_name[j]+'_93.5_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_W2,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_090-WMAP_W3_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmow3'+subtractmaps_name[j]+'_93.5_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_W3,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)
			smoothmap(indirectory,outdirectory,'CG_090-WMAP_W4_IQU_n0512_v1.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9cosmow4'+subtractmaps_name[j]+'_93.5_512_2023_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_W4,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True,usehealpixfits=True)

# EOF
