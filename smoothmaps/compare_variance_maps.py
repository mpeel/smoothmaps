from smoothmap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits
import os
import scipy.io as io
import matplotlib.pyplot as plt

indirectory = '/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0/'

namestrings = ['QUIJOTEMFI1_11.0_2021','QUIJOTEMFI1_13.0_2021','QUIJOTEMFI3_11.0_2021','QUIJOTEMFI3_13.0_2021','QUIJOTEMFI2_17.0_2021','QUIJOTEMFI2_19.0_2021','QUIJOTEMFI4_17.0_2021','QUIJOTEMFI4_19.0_2021']
for i in range(0,len(namestrings)):
	if i == 0 or i == 1:
		null = 0
	else:
		normal = hp.read_map(indirectory+'512_60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits',field=None)
		mask = np.ones(len(normal[0]))
		mask[normal[0] == hp.UNSEEN] = 0
		hp.mollview(mask)
		plt.savefig('qt_mask.png')
		plt.clf()
		sim = hp.read_map(indirectory+'512_60.0smoothed_'+namestrings[i].replace('_2021','_2021simnoise')+'_mKCMBunits.fits',field=None)
		normal[normal == 0.0] = hp.UNSEEN
		sim[sim == 0.0] = hp.UNSEEN
		hp.mollview(np.sqrt(normal[3])*mask,norm='hist')#,max=np.median(np.sqrt(normal[3])*10))
		plt.savefig(namestrings[i]+'_I_norm.png')
		plt.close()
		plt.clf()
		hp.mollview(np.sqrt(sim[3])*mask,norm='hist')#,max=np.median(np.sqrt(sim[3])*10))
		plt.savefig(namestrings[i]+'_I_sim.png')
		plt.clf()
		hp.mollview(sim[3]/normal[3]*mask)
		plt.savefig(namestrings[i]+'_I_ratio.png')
		plt.close()
		plt.clf()
		hp.mollview(np.sqrt(normal[4])*mask,norm='hist')#,max=np.median(np.sqrt(normal[4])*10))
		plt.savefig(namestrings[i]+'_Q_norm.png')
		plt.clf()
		hp.mollview(np.sqrt(sim[4])*mask,norm='hist')#,max=np.median(np.sqrt(sim[4])*10))
		plt.savefig(namestrings[i]+'_Q_sim.png')
		plt.clf()
		hp.mollview(sim[4]/normal[4]*mask)
		plt.savefig(namestrings[i]+'_Q_ratio.png')
		plt.close()
		plt.clf()
		print(namestrings[i])
		print(np.median(sim[3][mask==1]/normal[3][mask==1]))
		print(np.median(sim[4][mask==1]/normal[4][mask==1]))
		print(np.median(sim[6][mask==1]/normal[6][mask==1]))
		print(np.median(np.sqrt(sim[3][mask==1])/np.sqrt(normal[3][mask==1])))
		print(np.median(np.sqrt(sim[4][mask==1])/np.sqrt(normal[4][mask==1])))
		print(np.median(np.sqrt(sim[6][mask==1])/np.sqrt(normal[6][mask==1])))

		plt.hist(np.log10(sim[3][mask==1]/normal[3][mask==1]))
		plt.yscale('log')
		plt.savefig(namestrings[i]+'_I_hist.png')
		plt.close()
		plt.clf()
		plt.hist(np.log10(sim[4][mask==1]/normal[4][mask==1]))
		plt.yscale('log')
		plt.savefig(namestrings[i]+'_Q_hist.png')
		plt.close()
		plt.clf()
		plt.hist(np.log10(sim[6][mask==1]/normal[6][mask==1]))
		plt.yscale('log')
		plt.savefig(namestrings[i]+'_U_hist.png')
		plt.close()
		plt.clf()

		plt.hist(np.log10(normal[3][mask==1]))
		plt.hist(np.log10(sim[3][mask==1]))
		plt.yscale('log')
		plt.savefig(namestrings[i]+'_I_hist_comp.png')
		plt.close()
		plt.clf()
