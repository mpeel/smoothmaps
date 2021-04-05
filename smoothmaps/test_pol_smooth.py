import numpy as np
# import numba
import healpy as hp
from smoothmap import smoothmap, conv_nobs_variance_map
import astropy.io.fits as fits
import os
import scipy.io as io
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap

nside=512

inmap = '/Users/mpeel/Documents/maps/quijote_202103/reform/mfi_mar2021_11.0_3_smth.fits'
maps = hp.read_map(inmap,field=(0,1,2))
hp.gnomview(maps[1],reso=5,rot=[30,0])
plt.savefig('/Users/mpeel/Documents/maps/quijote_202103/smooth/q1.png')
plt.clf()
maps[0][maps[0][:] == hp.UNSEEN] = 0.0
maps[1][maps[1][:] == hp.UNSEEN] = 0.0
maps[2][maps[2][:] == hp.UNSEEN] = 0.0
inmap = '/Users/mpeel/Documents/maps/quijote_202103/smooth/512_60.0smoothed_QUIJOTEMFI3_11.0_2021_mKCMBunits.fits'
# inmap = '/Users/mpeel/Documents/maps/quijote_202103/smooth_old/512_60.0smoothed_QUIJOTEMFI3_11.0_2021_mKCMBunits.fits'
maps2 = hp.read_map(inmap,field=(0,1,2))
hp.gnomview(maps2[1],reso=5,rot=[30,0])
plt.savefig('/Users/mpeel/Documents/maps/quijote_202103/smooth/q.png')
plt.clf()
maps2[0][maps2[0][:] == hp.UNSEEN] = 0.0
maps2[1][maps2[1][:] == hp.UNSEEN] = 0.0
maps2[2][maps2[2][:] == hp.UNSEEN] = 0.0

hp.mollview((maps[1]-maps2[1])/maps[1],min=-1,max=1)
plt.savefig('/Users/mpeel/Documents/maps/quijote_202103/smooth/compare_q.png')
plt.clf()
hp.mollview((maps[2]-maps2[2])/maps[2],min=-1,max=1)
plt.savefig('/Users/mpeel/Documents/maps/quijote_202103/smooth/compare_u.png')
plt.clf()
exit()

# Window functions
wf = io.readsav('/Users/mpeel/Documents/maps/quijote_mfi/mfi_blconv_wl.sav')
# wfq = wf['wl_mfi'].T
wfq = wf['bl'].T
wfq_l = range(0,len(wfq[0][0]))

# maps = hp.read_map('test_smooth.fits',field=(0,1,2))
# hp.mollview(maps[1])
# plt.savefig('test_1.png')
# plt.clf()

# inmap = '/Users/mpeel/Documents/maps/quijote_202103/reform/mfi_mar2021_13.0_3.fits'
inmap = '/Users/mpeel/Documents/maps/quijote_202103/reform/mfi_mar2021_11.0_3_smth.fits'
# inmap = '/Users/mpeel/Documents/maps/wmap9/wmap_band_iqumap_r9_9yr_K_v5.fits'
maps = hp.read_map(inmap,field=(0,1,2))
# hp.mollview(maps[1])
# plt.savefig('test_2.png')
# plt.clf()
# exit()

maps = hp.read_map(inmap,field=(0,1,2))
cmap_colors = get_cmap('binary', 256)(np.linspace(0, 1, 256))
cmap_colors[..., 3] = 0.4  # Make colormap partially transparent
cmap = ListedColormap(cmap_colors)
lic = hp.line_integral_convolution(maps[1],maps[2])
lic = hp.smoothing(lic, np.deg2rad(0.5))
hp.mollview(np.sqrt(maps[1]**2 + maps[2]**2), cbar=True,unit='mK_CMB',max=1)
hp.mollview(lic, cmap=cmap,cbar=False, reuse_axes=True, title='MFI 11GHz')
plt.savefig('test_lic.png')
exit()

alms = hp.map2alm(maps)
for i in range(0,3):
	# alms[i] = hp.sphtfunc.smoothalm(alms[i], beam_window=wfq[2][1],pol=False)
	alms[i] = hp.sphtfunc.smoothalm(alms[i], fwhm=1.0*np.pi/180.0,pol=False)
newmap = hp.alm2map(alms, nside,verbose=False)

for i in range(0,3):
	# alms[i] = hp.sphtfunc.smoothalm(alms[i], beam_window=wfq[2][1],pol=False)
	alms[i] = hp.map2alm(maps[i])
	alms[i] = hp.sphtfunc.smoothalm(alms[i], fwhm=1.0*np.pi/180.0,pol=False)
newmap2 = hp.alm2map(alms, nside,verbose=False)

maps = hp.read_map(inmap,field=(0,1,2))
alms = hp.map2alm(maps,pol=False)
# alms = hp.sphtfunc.smoothalm(alms, beam_window=wfq[2][1],pol=False)
alms = hp.sphtfunc.smoothalm(alms, fwhm=1.0*np.pi/180.0,pol=False)
newmap2 = hp.alm2map(alms, nside,verbose=False,pol=False)

maps = hp.read_map(inmap,field=(0,1,2))
alms = hp.map2alm(maps,pol=True)
# alms = hp.sphtfunc.smoothalm(alms, beam_window=wfq[2][1],pol=True)
alms = hp.sphtfunc.smoothalm(alms, fwhm=1.0*np.pi/180.0,pol=True)
newmap3 = hp.alm2map(alms, nside,verbose=False)

maps = hp.read_map('test_smooth_wmap.fits',field=(0,1,2))
maps[newmap3==hp.UNSEEN]=hp.UNSEEN

hp.mollview(newmap[0]-newmap2[0])
plt.savefig('test_w_method12_I.png')
plt.clf()
hp.mollview(newmap3[0]-newmap2[0])
plt.savefig('test_w_method23_I.png')
plt.clf()
hp.mollview(newmap[0]-maps[0])
plt.savefig('test_w_method14_I.png')
plt.clf()
hp.mollview(newmap3[0]-maps[0])
plt.savefig('test_w_method34_I.png')
# plt.clf()

hp.mollview(newmap[1]-newmap2[1])
plt.savefig('test_w_method12_Q.png')
plt.clf()
hp.mollview(newmap3[1]-newmap2[1])
plt.savefig('test_w_method23_Q.png')
plt.clf()
hp.mollview(newmap[1]-maps[1])
plt.savefig('test_w_method14_Q.png')
plt.clf()
hp.mollview(newmap3[1]-maps[1])
plt.savefig('test_w_method34_Q.png')
plt.clf()

hp.mollview(newmap[2]-newmap2[2])
plt.savefig('test_w_method12_U.png')
plt.clf()
hp.mollview(newmap3[2]-newmap2[2])
plt.savefig('test_w_method23_U.png')
plt.clf()

hp.mollview(newmap[2]-maps[2])
plt.savefig('test_w_method14_U.png')
plt.clf()
hp.mollview(newmap3[2]-maps[2])
plt.savefig('test_w_method34_U.png')
plt.clf()
