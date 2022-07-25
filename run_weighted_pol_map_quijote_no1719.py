from weighted_pol_map import *
from fastcc.fastcc import *
import os

# General settings
nsides = [512]#[8, 16, 32, 64, 128, 256]
# indexes = [-3.0]
# indexes = [-2.9, -2.95, -3.0, -3.05, -3.1,-3.15, -3.2]
# indexes = [-1.8, -1.9, -2.0, -2.1]#[-2.2, -2.3, -2.4, -2.6]#, -2.5, -2.7, -2.8, -2.9, -3.0, -3.1, -3.2, -3.3, -3.4, -3.5]
indexes = [-3.1]
# Settings for this run file
doing_quijote = True

# Settings to pass to the code
use_halfrings = False
use_weights = False
use_reweight_by_rms = False
use_reweight_by_rms_method = 2 # 1 = ricardo, 2 = alberto
use_planck = False # This was for comparison only
use_cbass = False
freqs = [11.1,12.9]
normfreq = 10.0
use_planckwmap = True # Also combine with Planck/WMAP
# planckvers = ['2015','2015nobp','2018','2018nobp','2020']
use_extra_mask_for_qt = True
planckvers = ['2020']
only_wmap = False
only_planck = False
folder = 'weighted_wmap_planck_qu'
version='tqu_v1.5_noise_v1.0'
version_wmap = version + "_10k"
version_p18 = version + "_10k"
version_p20 = version + "_10k"
version_fdec = 'Fdec' # or set to '' to use the normal ones
doqu = True
minplots = True
maps_half1=[]
maps_half2=[]
separate_variance_maps=[]
statsmask = '/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_lowdec_nside512.fits'

for planckver in planckvers:
	for nside in nsides:
		for index in indexes:
			# QUIJOTE
			indirectory = '/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/'
			outdirectory = '/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_weighted_fdec/'
			date='202103'

			# # Set up QUIJOTE input
			# prefix='half1mfi'
			# maps_half1 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

			# prefix='half2mfi'
			# maps_half2 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

			usewei = True
			prefix='mfi'
			if usewei:
				maps = [str(nside)+'_60.0smoothed_QUIJOTEMFI3_11.0_2021_mKCMBunits.fits',str(nside)+'_60.0smoothed_QUIJOTEMFI3_13.0_2021_mKCMBunits.fits']
				prefix=str(nside)+'_60.0smoothed_quijotecombwei10_'+version+'_'+str(index)+'_no1719'
			else:
				maps = [str(nside)+'_60.0smoothed_QUIJOTEMFI3_11.0_2021simnoise_mKCMBunits.fits',str(nside)+'_60.0smoothed_QUIJOTEMFI3_13.0_2021simnoise_mKCMBunits.fits']
				prefix=str(nside)+'_60.0smoothed_quijotecomb10_'+version+'_'+str(index)+'_no1719'
			#_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',
			varianceindex=[[3,4,6,5], [3,4,6,5]]
			rescale_amp = np.ones(len(maps))
			rescale_variance = rescale_amp.copy()
			extra_mask = ''
			if use_extra_mask_for_qt == True:
				extra_mask = '/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_satband_nside512.fits'
			rescale_amp[0] *= fastcc('Q11',index+2.0,detector='Q311')
			rescale_amp[1] *= fastcc('Q13',index+2.0,detector='Q313')
			print(rescale_amp)
			if usewei:
				rescale_variance[0] *= 1.504
				rescale_variance[1] *= 1.401
			print(rescale_variance)


			if use_planckwmap:
				freqs = np.concatenate((freqs,[28.4, 44.1, 22.8, 33.0, 40.7]))
				varianceindex = np.concatenate((varianceindex,[[3,4,6,5], [3,4,6,5], [3,4,6,5], [3,4,6,5], [3,4,6,5]]))
				rescale_amp = np.concatenate((rescale_amp,np.ones(5)))
				rescale_variance = np.concatenate((rescale_variance, np.ones(5)))
				maps = maps + ['../planck2020_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_PlanckR4fullbeamnodp'+version_fdec+'Noise_28.4_1024_2020_mKCMBunits.fits','../planck2020_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_PlanckR4fullbeamnodp'+version_fdec+'Noise_44.1_1024_2020_mKCMBunits.fits','../wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beam'+version_fdec+'Noise_22.8_512_2013_mKCMBunits.fits','../wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beam'+version_fdec+'Noise_33.0_512_2013_mKCMBunits.fits','../wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beam'+version_fdec+'Noise_40.7_512_2013_mKCMBunits.fits']
				# print(maps)
				# exit()
			# 	# indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/'
			# 	# outdirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_weight/'

			# 	# Planck and WMAP colour corrections
				rescale_amp[2] *= fastcc('P30',index+2.0)
				rescale_amp[3] *= fastcc('P44',index+2.0)
				rescale_amp[4] *= fastcc('WK',index+2.0)
				rescale_amp[5] *= fastcc('WKa',index+2.0)
				rescale_amp[6] *= fastcc('WQ',index+2.0)
				prefix = prefix + '_wmapplanck'

			# This is only Planck+WMAP
			if False:
				freqs = [28.4, 44.1, 22.8, 33.0, 40.7]
				varianceindex = [[3,4,6,5], [3,4,6,5], [3,4,6,5], [3,4,6,5], [3,4,6,5]]
				rescale_amp = np.ones(5)
				rescale_variance = np.ones(5)
				maps = ['../planck2020_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_PlanckR4fullbeamnodp'+version_fdec+'Noise_28.4_1024_2020_mKCMBunits.fits','../planck2020_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_PlanckR4fullbeamnodp'+version_fdec+'Noise_44.1_1024_2020_mKCMBunits.fits','../wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beam'+version_fdec+'Noise_22.8_512_2013_mKCMBunits.fits','../wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beam'+version_fdec+'Noise_33.0_512_2013_mKCMBunits.fits','../wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beam'+version_fdec+'Noise_40.7_512_2013_mKCMBunits.fits']
				# print(maps)
				# exit()
			# 	# indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/'
			# 	# outdirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_weight/'

			# 	# Planck and WMAP colour corrections
				rescale_amp[0] *= fastcc('P30',index+2.0)
				rescale_amp[1] *= fastcc('P44',index+2.0)
				rescale_amp[2] *= fastcc('WK',index+2.0)
				rescale_amp[3] *= fastcc('WKa',index+2.0)
				rescale_amp[4] *= fastcc('WQ',index+2.0)
				prefix = 'wmapplanck'
			# 	apply_extra_mask=np.zeros(len(freqs))
			# 	apply_extra_mask[startindex+0] = 1
			# 	apply_extra_mask[startindex+1] = 1
			# 	extra_mask='compare_wmap_planck_polmap_mask.fits'

			apply_extra_mask=np.zeros(len(freqs))
			# Apply the extra mask to QUIJOTE?
			if use_extra_mask_for_qt == True:
				apply_extra_mask[0] = 1
				apply_extra_mask[1] = 1
				prefix = prefix + '_qtmask'

			# Make sure the output directory exists
			os.makedirs(outdirectory, exist_ok=True)
			# ... and run the weighted map!
			print(maps)
			print(prefix)
			weighted_pol_map(nside=nside,indirectory=indirectory,outdirectory=outdirectory,date=date,prefix=prefix,index=index,freqs=freqs,maps=maps,maps_half1=maps_half1,maps_half2=maps_half2,use_halfrings=use_halfrings,use_weights=use_weights,use_reweight_by_rms=use_reweight_by_rms,use_reweight_by_rms_method=use_reweight_by_rms_method,use_planck=use_planck,use_cbass=use_cbass,normfreq=normfreq,rescale_amp=rescale_amp,rescale_variance=rescale_variance,apply_extra_mask=apply_extra_mask,extra_mask=extra_mask,varianceindex=varianceindex,threshold=1.0,separate_variance_maps=separate_variance_maps,doqu=doqu,minplots=minplots,dodiffs=False,statsmask=statsmask)
