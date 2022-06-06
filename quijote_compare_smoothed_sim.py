import astropy.io.fits as fits
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.io as io

from smoothmap import *

# Window functions
directory = '/Users/mpeel/Documents/maps/quijote_202103/reform/'
wf = io.readsav(directory+'mfi_blconv_wl.sav')
# wfq = wf['wl_mfi'].T
wfq = wf['bl'].T
wfq_l = range(0,len(wfq[0][0]))

directory = '/Users/mpeel/Desktop/'
infile = directory+'quijote_11GHz_horn3_0444_sm1deg.fits'

# inputfits = fits.open(infile)
# print(inputfits.info())
# cols = inputfits[1].columns
# col_names = cols.names
# nmaps = len(cols)
# maps = []
# for i in range(0,nmaps):
# 	maps.append(inputfits[1].data.field(i))

# alms = hp.map2alm(maps[0:3],pol=True)

inputfits = fits.open(infile)
print(inputfits.info())
print(len(inputfits[1].data.field(0)))
nside=512
alms = hp.map2alm([inputfits[1].data.field(0),inputfits[1].data.field(1),inputfits[1].data.field(2)],pol=True)
alms = hp.sphtfunc.smoothalm(alms, beam_window=wfq[2][0],pol=True)
newmap = hp.alm2map(alms, nside,verbose=False,pol=True)
hp.write_map(directory+'testmap.fits',newmap)
hp.mollview(inputfits[1].data.field(3) - newmap[0])
plt.savefig(directory+'compare_I2.png')
plt.clf()
hp.mollview(inputfits[1].data.field(4) - newmap[1])
plt.savefig(directory+'compare_Q2.png')
plt.clf()

exit()

# smoothmap(directory,directory,'quijote_11GHz_horn3_0444_sm1deg.fits','quijote_11GHz_horn3_0444_sm1deg_checksmth2.fits', 0.0,windowfunction=wfq[2][0],nside_out=512,units_in='mKCMB',units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=False)
# exit()

inputfits = fits.open(infile)
nmaps = len(inputfits[1].columns)
print(inputfits[0].header)
# exit()
maps = []
for i in range(0,nmaps):
	maps.append(inputfits[1].data.field(i))
maps[3][maps[3][:] == hp.UNSEEN] = 0.0
maps[4][maps[4][:] == hp.UNSEEN] = 0.0
maps[5][maps[5][:] == hp.UNSEEN] = 0.0
hp.mollview(maps[4])
plt.savefig(directory+'test_Q.png')
plt.clf()
hp.mollview(maps[0])
plt.savefig(directory+'test_I_raw.png')
plt.clf()
hp.mollview(maps[3])
plt.savefig(directory+'test_I_mateo.png')
plt.clf()
hp.mollview(maps[3]-np.median(maps[3]))
plt.savefig(directory+'test_I_mateo_medsub.png')
plt.clf()
infile2 = directory+'quijote_11GHz_horn3_0444_sm1deg_checksmth2.fits'
inputfits2 = fits.open(infile2)
nmaps = len(inputfits[1].columns)
maps2 = []
for i in range(0,nmaps):
	maps2.append(inputfits2[1].data.field(i))
maps2[0][maps2[0][:] == hp.UNSEEN] = 0.0
maps2[1][maps2[1][:] == hp.UNSEEN] = 0.0
maps2[2][maps2[2][:] == hp.UNSEEN] = 0.0
maps2[0][maps2[0][:] < -1e10] = 0.0
maps2[1][maps2[1][:] < -1e10] = 0.0
maps2[2][maps2[2][:] < -1e10] = 0.0
hp.mollview(maps2[1])
plt.savefig(directory+'test_Q2.png')
plt.clf()
hp.mollview(maps2[0])
plt.savefig(directory+'test_I_mike.png')
plt.clf()
hp.mollview(maps2[0]-np.median(maps2[0]))
plt.savefig(directory+'test_I_mike_medsub.png')
plt.clf()

# maps = hp.read_map(infile,field=None)
# maps2 = hp.read_map(infile2,fields=(1,2,3))


hp.mollview(maps[3] - maps2[0],title='I')
plt.savefig(directory+'compare_I.png')
plt.clf()
hp.mollview(maps[4] - maps2[1],title='Q')
plt.savefig(directory+'compare_Q.png')
plt.clf()
hp.mollview(maps[4] - maps2[1],title='Q',norm='hist')
plt.savefig(directory+'compare_Q_hist.png')
plt.clf()
hp.mollview(maps[5] - maps2[2],title='U')
plt.savefig(directory+'compare_U.png')
plt.clf()
