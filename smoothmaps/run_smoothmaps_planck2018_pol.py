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
outdirectory = directory+"planck2018_tqu_v1.5/"
os.makedirs(outdirectory, exist_ok=True)

beamtf_p30 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',28)
beamtf_p44 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',29)
beamtf_p70 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',30)

beamtf_p100 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_100x100.fits')
beamtf_p143 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_143x143.fits')
beamtf_p217 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_217x217.fits')
beamtf_p353 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_353x353.fits')
beamtf_p545 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_545x545.fits')
beamtf_p857 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_857x857.fits')

numnside = len(output_nside)
for i in range(0,numnside):
	# Smooth CMB maps - just saved, not used below at all.
	smoothmap(directory+'planck2018/',outdirectory,'COM_CMB_IQU-commander_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBCommander_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB',do_pol_combined=True)
	smoothmap(directory+'planck2018/',outdirectory,'COM_CMB_IQU-nilc_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBNILC_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB',do_pol_combined=True)
	smoothmap(directory+'planck2018/',outdirectory,'COM_CMB_IQU-sevem_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBSevem_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB',do_pol_combined=True)
	smoothmap(directory+'planck2018/',outdirectory,'COM_CMB_IQU-smica_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBSmica_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB',do_pol_combined=True)
	smoothmap(directory+'planck2018/',outdirectory,'COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBSmica-nosz_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB',do_pol_combined=True)

	subtractmaps = ['']
	subtractmaps_name = ['']
	numsubtract = len(subtractmaps)
	for j in range(0,numsubtract):
		if output_nside[i] <= 1024:
			# LFI nominal
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_030_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeamnobp'+subtractmaps_name[j]+'_28.4_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p30,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_044_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeamnobp'+subtractmaps_name[j]+'_44.1_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p44,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_070_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeamnobp'+subtractmaps_name[j]+'_70.4_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p70,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

			# LFI bandpass corrected
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_28.4_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p30,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_44.1_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p44,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_70.4_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p70,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

			# LFI deconvolved
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_030-deconvolved-040-IQUSS_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3dec40'+subtractmaps_name[j]+'_28.4_1024_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-40.0**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_044-deconvolved-030-IQUSS_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3dec30'+subtractmaps_name[j]+'_44.1_1024_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-30.0**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
			smoothmap(directory+'planck2018/',outdirectory,'LFI_SkyMap_070-deconvolved-020-IQUSS_1024_R3.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3dec20'+subtractmaps_name[j]+'_70.4_1024_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-20.0**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

		# HFI polarisation channels
		smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_100_2048_R3.01_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_100_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p100,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
		smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_143_2048_R3.01_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_143_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p143,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
		smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_217_2048_R3.01_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_217_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p217,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)
		smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_353-psb_2048_R3.01_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_353_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p353,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,do_pol_combined=True)

		# HFI unpolarised channels
		if subtractmaps_name[j] == '':
			smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_545_2048_R3.01_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_545_2048_2018_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p545,smoothvariance=smoothvariance,do_pol_combined=False)
			smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_857_2048_R3.01_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_857_2048_2018_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p857,smoothvariance=smoothvariance,do_pol_combined=False)


# EOF