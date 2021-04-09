from smoothmap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits
import os
import scipy.io as io
import matplotlib.pyplot as plt

# Window functions
directory = '/Users/mpeel/Documents/maps/quijote_202103/reform/'
wf = io.readsav(directory+'mfi_blconv_wl.sav')
# wfq = wf['wl_mfi'].T
wfq = wf['bl'].T
wfq_l = range(0,len(wfq[0][0]))

directory = '/Users/mpeel/Desktop/'
infile = directory+'quijote_11GHz_horn3_0001_sm1deg.fits'
inputfits = fits.open(infile)
nmaps = len(inputfits[1].columns)
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
infile2 = directory+'quijote_11GHz_horn3_0001_sm1deg_checksmth.fits'
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


hp.mollview(maps[3] - maps2[0])
plt.savefig(directory+'compare_I.png')
plt.clf()
hp.mollview(maps[4] - maps2[1])
plt.savefig(directory+'compare_Q.png')
plt.clf()
hp.mollview(maps[5] - maps2[2])
plt.savefig(directory+'compare_U.png')
plt.clf()


# smoothmap(directory,directory,'quijote_11GHz_horn3_0001_sm1deg.fits','quijote_11GHz_horn3_0001_sm1deg_checksmth.fits', 60.0,windowfunction=wfq[2][0],nside_out=512,units_out='mKCMB',use_precomputed_wf=True,do_pol_combined=True,useunseen=True,)

