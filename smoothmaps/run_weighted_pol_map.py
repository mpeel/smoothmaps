from weighted_pol_map import *
from astrocode.colourcorrections.fastcc import *
import os

# General settings
nsides = [8, 16, 32, 64, 128, 256]
indexes = [-2.9, -2.95, -3.0, -3.05, -3.1]

# Settings for this run file
doing_quijote = False

# Settings to pass to the code
use_halfrings = False
use_weights = False
use_reweight_by_rms = False
use_reweight_by_rms_method = 2 # 1 = ricardo, 2 = alberto
use_planck = False # This was for comparison only
use_cbass = False
freqs = [17,19,11,13,17,19]
# normfreq = 10.0
normfreq = 28.4
# planckver = '2015'
# planckver = '2015nobp'
# planckver = '2018'
# planckver = '2018nobp'
# planckver = '2018dec'
# planckver = '2020'
planckvers = ['2015','2015nobp','2018','2018nobp','2020']
only_wmap = False
only_planck = False
version='tqu_v1.4_noise_v0.8'

for planckver in planckvers:
	for nside in nsides:
		for index in indexes:
			if doing_quijote == False:
				date = ''
				maps_half1=[]
				maps_half2=[]
				separate_variance_maps=[]

				freqs = [28.4, 44.1, 22.8, 33.0, 40.7]
				rescale_amp = np.ones(len(freqs))

				indirectory = '/Users/mpeel/Documents/maps/'
				if planckver == '2015':
					outdirectory = indirectory+'weighted_wmap_planck_'+version+'/'
					prefix=str(nside)+'_60.0smoothed_wmap9_planck2015_'+version+'_'+str(index)
					maps = ['planck2015_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR2fullbeambpcorrNoise_28.4_256_2015_mKCMBunits.fits','planck2015_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR2fullbeambpcorrNoise_44.1_256_2015_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
					varianceindex=[[3,4,5], [3,4,5], [3,4,5], [3,4,5], [3,4,5]]
				elif planckver == '2015nobp':
					outdirectory = indirectory+'weighted_wmap_planck_'+version+'/'
					prefix=str(nside)+'_60.0smoothed_wmap9_planck2015nobp_'+version+'_'+str(index)
					maps = ['planck2015_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR2fullbeamNoise_28.4_1024_2015_mKCMBunits.fits','planck2015_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR2fullbeamNoise_44.1_1024_2015_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
					varianceindex=[[3,4,5], [3,4,5], [3,4,5], [3,4,5], [3,4,5]]
				elif planckver == '2018':
					outdirectory = indirectory+'weighted_wmap_planck_'+version+'/'
					prefix=str(nside)+'_60.0smoothed_wmap9_planck2018_'+version+'_'+str(index)
					maps = ['planck2018_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR3fullbeamNoise_28.4_1024_2018_mKCMBunits.fits','planck2018_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR3fullbeamNoise_44.1_1024_2018_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
					varianceindex=[[3,4,5], [3,4,5], [3,4,5], [3,4,5], [3,4,5]]
				elif planckver == '2018nobp':
					outdirectory = indirectory+'weighted_wmap_planck_'+version+'/'
					prefix=str(nside)+'_60.0smoothed_wmap9_planck2018nobp_'+version+'_'+str(index)
					maps = ['planck2018_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR3fullbeamnobpNoise_28.4_1024_2018_mKCMBunits.fits','planck2018_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR3fullbeamnobpNoise_44.1_1024_2018_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
					varianceindex=[[3,4,5], [3,4,5], [3,4,5], [3,4,5], [3,4,5]]
				elif planckver == '2018dec':
					outdirectory = indirectory+'weighted_wmap_planck_'+version+'/'
					prefix=str(nside)+'_60.0smoothed_wmap9_planck2018dec_'+version+'_'+str(index)
					maps = ['planck2018_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR3decNoise_dec20_70.4_1024_2018_mKCMBunits.fits','planck2018_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR3decNoise_dec30_44.1_1024_2018_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
					varianceindex=[[3,4,5], [3,4,5], [3,4,5], [3,4,5], [3,4,5]]
				elif planckver == '2020':
					outdirectory = indirectory+'weighted_wmap_planck_'+version+'/'
					prefix=str(nside)+'_60.0smoothed_wmap9_planck2020_'+version+'_'+str(index)
					maps = ['planck2020_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits','planck2020_'+version+'/'+str(nside)+'_60.0smoothed_PlanckR4fullbeamnodpNoise_44.1_1024_2020_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
					varianceindex=[[3,4,5], [3,4,5], [3,4,5], [3,4,5], [3,4,5]]

				rescale_variance = rescale_amp.copy()

				# prefix='wmap9_planck2015_tqu_10ghz'
				# freqs = [28.4, 44.1, 22.8, 33.0, 40.7]
				# rescale_amp = [1.0, 1.0, 1.0, 1.0, 1.0]
				# maps = ['512_60.0smoothed_PlanckR2fullbeambpcorr_28.4_256_2015_mKCMBunits.fits','512_60.0smoothed_PlanckR2fullbeambpcorr_44.1_256_2015_mKCMBunits.fits','../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']
				# indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2015_tqu/'
				# outdirectory = '/Users/mpeel/Documents/maps/wmap9_planck2015_weight/'

				# Planck and WMAP colour corrections
				rescale_amp[0] *= fastcc('30',index+2.0)
				rescale_amp[1] *= fastcc('44',index+2.0)
				rescale_amp[2] *= fastcc('K',index+2.0)
				rescale_amp[3] *= fastcc('Ka',index+2.0)
				rescale_amp[4] *= fastcc('Q',index+2.0)
				print(rescale_amp)

				# Rescale the WMAP variances as the wrong sigma_0 was used to generate them.
				# No longer needed for noise realisations of v0.8
				# rescale_variance[2] *= (1.435/1.429)**2
				# rescale_variance[3] *= (1.472/1.466)**2
				# rescale_variance[4] *= (2.197/2.188)**2
				# print(rescale_variance)
				# exit()

				# Test rescaling the variances by the difference in beam size
				# rescale_variance[0] *= (32.0/60.0)**4
				# rescale_variance[1] *= (27.0/60.0)**4
				# rescale_variance[2] *= (49.0/60.0)**4
				# rescale_variance[3] *= (40.0/60.0)**4
				# rescale_variance[4] *= (31.0/60.0)**4

				print(rescale_variance)

				apply_extra_mask=np.zeros(len(freqs))
				extra_mask=''
				# apply_extra_mask[0] = 1
				# apply_extra_mask[1] = 1
				# extra_mask='compare_wmap_planck_polmap_mask.fits'

				if only_wmap:
					prefix='wmap9_'+version
					freqs = [22.8, 33.0, 40.7]
					maps = ['wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','wmap9_'+version+'/'+str(nside)+'_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
					rescale_amp = np.ones(len(freqs))
					rescale_variance = np.ones(len(freqs))
					rescale_amp[0] *= fastcc('K',index+2.0)
					rescale_amp[1] *= fastcc('Ka',index+2.0)
					rescale_amp[2] *= fastcc('Q',index+2.0)
					varianceindex=[[3,4,5], [3,4,5], [3,4,5]]

				if only_planck:
					maps = maps[0:2]
					prefix = prefix.replace('wmap9_','')
				# prefix='planck2018'
				# freqs = [28.4, 44.1]
				# maps = ['512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits']

				# prefix='wmap9_planck2018_10ghz'
				# normfreq = 10.0
				# indirectory = '/Users/mpeel/Documents/maps/quijote_201907/smooth/'
				# outdirectory = '/Users/mpeel/Documents/maps/quijote_201907/weighted/'

				# maps = ['../../wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']
			else:
				# QUIJOTE
				indirectory = '/Users/mpeel/Documents/maps/quijote_201907/smooth/'
				outdirectory = '/Users/mpeel/Documents/maps/quijote_201907/weighted/'
				date='201907'

				# indirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
				# outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/analyse/'
				# date='201905'

				# Set up QUIJOTE input
				prefix='half1mfi'
				maps_half1 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

				prefix='half2mfi'
				maps_half2 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

				prefix='mfi'
				maps = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

				prefix='mfi_10ghz'
				rescale_amp = np.ones(len(maps))
				rescale_variance = rescale_amp.copy()
				apply_extra_mask=np.zeros(len(freqs))
				extra_mask = ''


				if use_planck:
					startindex = len(freqs)
					freqs = [17,19,11,13,17,19,28.4, 44.1, 22.8, 33.0, 40.7]
					rescale_amp = np.ones(len(freqs))
					rescale_variance = np.ones(len(freqs))
					prefix='mfi'
					maps = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','../../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']
					print(maps)
					# indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/'
					# outdirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_weight/'

					prefix='mfi_planck_10ghz'
				

					# Planck and WMAP colour corrections
					rescale_amp[startindex+0] *= fastcc('30',index+2.0)
					rescale_amp[startindex+1] *= fastcc('44',index+2.0)
					rescale_amp[startindex+2] *= fastcc('K',index+2.0)
					rescale_amp[startindex+3] *= fastcc('Ka',index+2.0)
					rescale_amp[startindex+4] *= fastcc('Q',index+2.0)
					print(rescale_amp)

					# Rescale the WMAP variances as the wrong sigma_0 was used to generate them.
					# No longer needed for noise realisations of v0.8
					# rescale_variance[startindex+2] *= (1.435/1.429)**2
					# rescale_variance[startindex+3] *= (1.472/1.466)**2
					# rescale_variance[startindex+4] *= (2.197/2.188)**2
					# print(rescale_variance)
					# exit()

					apply_extra_mask=np.zeros(len(freqs))
					apply_extra_mask[startindex+0] = 1
					apply_extra_mask[startindex+1] = 1
					extra_mask='compare_wmap_planck_polmap_mask.fits'

			# Make sure the output directory exists
			os.makedirs(outdirectory, exist_ok=True)
			# ... and run the weighted map!
			print(maps)
			weighted_pol_map(nside=nside,indirectory=indirectory,outdirectory=outdirectory,date=date,prefix=prefix,index=index,freqs=freqs,maps=maps,maps_half1=maps_half1,maps_half2=maps_half2,use_halfrings=use_halfrings,use_weights=use_weights,use_reweight_by_rms=use_reweight_by_rms,use_reweight_by_rms_method=use_reweight_by_rms_method,use_planck=use_planck,use_cbass=use_cbass,normfreq=normfreq,rescale_amp=rescale_amp,rescale_variance=rescale_variance,apply_extra_mask=apply_extra_mask,extra_mask=extra_mask,varianceindex=varianceindex,threshold=0.1,separate_variance_maps=separate_variance_maps)
