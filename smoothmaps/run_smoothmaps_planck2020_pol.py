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

output_resolutions = [60.0, 20.0]
output_nside = np.asarray([2048, 1024, 512, 256, 128, 64, 32, 16, 8])
smoothvariance = False

# directory = '/Users/mpeel/Documents/maps/'
# directory = '/scratch1/mpeel/maps/'
directory = '/share/nas_cbassarc/mpeel/'
outdirectory = directory+"planck2020_tqu_v1.5/"
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

dipole = 'dipole_nside2048.fits'
subtractmap_units='KCMB'

for output_resolution in output_resolutions:
	numnside = len(output_nside)
	for i in range(0,numnside):
		subtractmaps = [dipole,'']
		subtractmaps_name = ['nodp','']
		numsubtract = len(subtractmaps)
		for j in range(0,numsubtract):

			if output_resolution > 30.0:
				smoothmap(directory+'planck2020/',outdirectory,'LFI_SkyMap_030_1024_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_28.4_1024_2020_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p30,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units,do_pol_combined=True)
				smoothmap(directory+'planck2020/',outdirectory,'LFI_SkyMap_044_1024_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_44.1_1024_2020_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p44,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units,do_pol_combined=True)
				smoothmap(directory+'planck2020/',outdirectory,'LFI_SkyMap_070_1024_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_70.4_1024_2020_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p70,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units,do_pol_combined=True)

			smoothmap(directory+'planck2020/',outdirectory,'HFI_SkyMap_100_2048_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_100_2048_2020_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p100,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units,do_pol_combined=True)

			smoothmap(directory+'planck2020/',outdirectory,'HFI_SkyMap_143_2048_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_143_2048_2020_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p143,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units,do_pol_combined=True)

			smoothmap(directory+'planck2020/',outdirectory,'HFI_SkyMap_217_2048_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_217_2048_2020_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p217,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units,do_pol_combined=True)

			smoothmap(directory+'planck2020/',outdirectory,'HFI_SkyMap_353_2048_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_353_2048_2020_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p353,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units,do_pol_combined=True)

			if subtractmaps_name[j] == '':
				smoothmap(directory+'planck2020/',outdirectory,'HFI_SkyMap_545_2048_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_545_2048_2020_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p545,smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units)

				smoothmap(directory+'planck2020/',outdirectory,'HFI_SkyMap_857_2048_R4.00_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR4fullbeam'+subtractmaps_name[j]+'_857_2048_2020_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p857,smoothvariance=smoothvariance,usehealpixfits=True,subtractmap_units=subtractmap_units)

# EOF