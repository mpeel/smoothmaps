import astropy.io.fits as fits
import healpy as hp
import numpy as np
import os

from smoothmap import *

def get_hfi_beam(FITSfile):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print(header)
	print(data.names) # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

output_resolution = 60.0
output_nside = [2048, 1024, 512, 256, 128, 64, 32, 16, 8]
smoothvariance = False

# directory = '/Users/mpeel/Documents/maps/'
# directory = '/scratch1/mpeel/maps/'
directory = '/share/nas_cbassarc/mpeel/'
outdirectory = directory+"wmap9_tqu_v1.5/"
os.makedirs(outdirectory, exist_ok=True)

beamtf_K = np.loadtxt(directory+'wmap9/wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_V = np.loadtxt(directory+'wmap9/wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_W = np.loadtxt(directory+'wmap9/wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

numnside = len(output_nside)
for i in range(0,numnside):
	subtractmaps = ['']
	subtractmaps_name = ['']
	numsubtract = len(subtractmaps)
	for j in range(0,numsubtract):

		if output_nside[i] <= 512:
			# NB: these only use sigma_0 and Nobs for the intensity maps!
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_K_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_22.8_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK',windowfunction=beamtf_K,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_33.0_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK',windowfunction=beamtf_Ka,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_40.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK',windowfunction=beamtf_Q,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_V_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_60.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK',windowfunction=beamtf_V,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_W_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_93.5_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK',windowfunction=beamtf_W,subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

# EOF
