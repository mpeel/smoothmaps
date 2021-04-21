from smoothmap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits
import os

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
outdirectory = directory+"planck2015_tqu_v1.5/"
os.makedirs(outdirectory, exist_ok=True)

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

numnside = len(output_nside)
for i in range(0,numnside):
	subtractmaps = ['']
	subtractmaps_name = ['']
	numsubtract = len(subtractmaps)
	for j in range(0,numsubtract):
		# NB: These are 1° smoothed already!!
		if output_nside[i] <= 256:
			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_030-BPassCorrected_0256_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeambpcorr'+subtractmaps_name[j]+'_28.4_256_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-(60.0)**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_044-BPassCorrected_0256_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeambpcorr'+subtractmaps_name[j]+'_44.1_256_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-(60.0)**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_070-BPassCorrected_0256_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeambpcorr'+subtractmaps_name[j]+'_70.4_256_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-(60.0)**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

		if output_nside[i] <= 1024:
			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_030_1024_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_28.4_1024_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p30,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_044_1024_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_44.1_1024_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p44,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
		
		smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_070_2048_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_70.4_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p70,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_100_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_100_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p100,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_143_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_143_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p143,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_217_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_217_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p217,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_353_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_353_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p353,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

		# Only want these for the non-CMB subtracted sets.
		if subtractmaps_name[j] == '':
			smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_545_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_545_2048_2015_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p545,smoothvariance=smoothvariance,do_pol_combined=False)
			smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_857_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_857_2048_2015_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p857,smoothvariance=smoothvariance,do_pol_combined=False)

# EOF