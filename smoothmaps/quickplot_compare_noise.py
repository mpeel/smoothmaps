import matplotlib.pyplot as plt
import healpy as hp
import numpy as np

# Nside=512
map1_info = {'filename': '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits', 'name': 'WMAP33','rescale':1.0}
map2_info = {'filename': '/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits','name': 'Planck2020_28','rescale':1.0}
outputname = 'quickplot_compare_noise_wmapka_planck30.png'
range=[0,10.0]

# Nside=64
map1_info = {'filename': '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits', 'name': 'WMAP33','rescale':1.0}
map2_info = {'filename': '/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits','name': 'Planck2020_28','rescale':1.0}
outputname = 'quickplot_compare_noise_wmapka_planck30_64.png'
range=[0,60.0]

# To K-band
# map1_info = {'filename': '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits', 'name': 'WMAP22','rescale':1.0}
# map2_info = {'filename': '/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits','name': 'Planck2020_28','rescale':1.0}
# outputname = 'quickplot_compare_noise_wmapk_planck30_64.png'
# range=[0,0.05]

# map1_info = {'filename': '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits', 'name': 'WMAP33','rescale':1.0}
# map2_info = {'filename': '/Users/mpeel/Documents/maps/planck2020/LFI_SkyMap_030_1024_R4.00_full.fits','name': 'Planck2020_28','rescale':1e3}
# outputname = 'quickplot_compare_noise_wmapka_planck30_nosmooth.png'
# range=[0,0.1]


nside = 512

map1 = hp.read_map(map1_info['filename'],field=None)
map1 = hp.ud_grade(map1,nside_out=nside,power=2)
map2 = hp.read_map(map2_info['filename'],field=None)
map2 = hp.ud_grade(map2,nside_out=nside,power=2)
map1[4] *= (map1_info['rescale'])**2
map2[4] *= (map2_info['rescale'])**2
map1_med = np.sqrt(np.median(map1[4]))*1e3
print(map1_med)
map2_med = np.sqrt(np.median(map2[4]))*1e3
print(map2_med)
map1_mean = np.sqrt(np.mean(map1[4]))*1e3
map2_mean = np.sqrt(np.mean(map2[4]))*1e3
print(map1_med/map2_med)
print(map1_mean/map2_mean)

float_formatter = "{:.4f}".format

map1_n, map1_bins, null = plt.hist(np.sqrt(map1[4])*1e3,50, facecolor='blue', alpha=0.5,range=range,label=map1_info['name'] + ' med/mean ' + float_formatter(map1_med) + '/'+float_formatter(map1_mean))
map2_n, map2_bins, null = plt.hist(np.sqrt(map2[4])*1e3,50, facecolor='orange', alpha=0.5,range=range,label=map2_info['name'] + ' med/mean ' + float_formatter(map2_med) + '/'+float_formatter(map2_mean))
plt.xlabel('Noise level (uK_CMB)')
plt.ylabel('Counts')
plt.legend()
plt.savefig(outputname)
