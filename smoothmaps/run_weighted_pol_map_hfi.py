from weighted_pol_map import *
from astrocode.colourcorrections.fastcc import *
from astrocode.fitspectrum.spectra import *
# General settings
nside = 2048
index = 0.0
beta_d = 1.53 # From Planck 2018 IV
T_d = 19.6 # Ditto

# Settings for this run file
doing_quijote = False
# Settings to pass to the code
use_halfrings = False
use_weights = False
use_reweight_by_rms = True
use_reweight_by_rms_method = 2 # 1 = ricardo, 2 = alberto
use_planck = False
use_cbass = False
# freqs = [17,19,11,13,17,19]#11,13,
# normfreq = 10.0
normfreq = 353.0

date = ''
maps_half1=[]
maps_half2=[]

# prefix='planck2020_hfi_tqu_60arcmin_v1'
prefix='planck2020_hfi_tqu_20arcmin_v1'
freqs = [353.0,217.0,143.0]
rescale_amp = [1.0, 1.0,1.0]#, 1.0, 1.0, 1.0]
rescale_variance = rescale_amp.copy()
const = get_spectrum_constants()
val_353 = thermaldust_comm(const, 353.0, 1.0, beta_d, T_d)*planckcorr(const, 353.0)
val_217 = thermaldust_comm(const, 217.0, 1.0, beta_d, T_d)*planckcorr(const, 217.0)
val_143 = thermaldust_comm(const, 143.0, 1.0, beta_d, T_d)*planckcorr(const, 143.0)
print(str(val_353))
print(str(val_353/val_217) + ' - ' + str((normfreq/217)**beta_d))
print(str(val_353/val_143) + ' - ' + str((normfreq/143)**beta_d))
# exit()
rescale_amp[1] = val_353/val_217
rescale_variance[1] = rescale_amp[1]**2
rescale_amp[2] = val_353/val_143
rescale_variance[2] = rescale_amp[2]**2
maps = ['2048_20.0smoothed_PlanckR4fullbeamNoise_353_2048_2020_mKCMBunits.fits','2048_20.0smoothed_PlanckR4fullbeamNoise_217_2048_2020_mKCMBunits.fits','2048_20.0smoothed_PlanckR4fullbeamNoise_143_2048_2020_mKCMBunits.fits']
# maps = ['2048_60.0smoothed_PlanckR4fullbeamNoise_353_2048_2020_mKCMBunits.fits','2048_60.0smoothed_PlanckR4fullbeamNoise_217_2048_2020_mKCMBunits.fits','2048_60.0smoothed_PlanckR4fullbeamNoise_143_2048_2020_mKCMBunits.fits']
indirectory = '/Users/mpeel/Documents/maps/planck2020_tqu_v1.4_noise_v0.8_20arcmin/'
# indirectory = '/Users/mpeel/Documents/maps/planck2020_tqu_v1.4_noise_v0.8/'
outdirectory = '/Users/mpeel/Documents/maps/weighted_wmap_planck_tqu_v1.4_noise_v0.8/'
varianceindex=[[3,4,5], [3,4,5], [3,4,5]]
apply_extra_mask=np.zeros(len(freqs))
apply_extra_mask[0] = 0
apply_extra_mask[1] = 0
extra_mask=''#compare_wmap_planck_polmap_mask.fits'
threshold=0.1

weighted_pol_map(nside=nside,indirectory=indirectory,outdirectory=outdirectory,date=date,prefix=prefix,index=index,freqs=freqs,maps=maps,maps_half1=maps_half1,maps_half2=maps_half2,use_halfrings=use_halfrings,use_weights=use_weights,use_reweight_by_rms=use_reweight_by_rms,use_reweight_by_rms_method=use_reweight_by_rms_method,use_planck=use_planck,use_cbass=use_cbass,normfreq=normfreq,rescale_amp=rescale_amp,rescale_variance=rescale_variance,apply_extra_mask=apply_extra_mask,extra_mask=extra_mask,threshold=threshold,varianceindex=varianceindex)
