#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick analysis of the MFI maps
# 
# Version history:
#
# 31-May-2019  M. Peel       Started
# 05-Jun-2019  M. Peel       Generalised to cope with multiple runs
# 21-Dec-2020  M. Peel		 Add QU support

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astrocode.colourcorrections.fastcc import *
import astropy.io.fits as fits

def noiserealisation(inputmap, numpixels):
    newmap = np.zeros(numpixels)
    newmap = np.random.normal(scale=1.0, size=numpixels) * inputmap
    return newmap

def weighted_pol_map(nside=512,indirectory='',outdirectory='',date='',prefix='',index=-3.0,freqs=[],maps=[],maps_half1=[],maps_half2=[],use_halfrings=False,use_weights=False,use_reweight_by_rms=False,use_reweight_by_rms_method=2,use_planck=True,use_cbass=False,normfreq=10.0,rescale_amp=[],rescale_variance=[],apply_extra_mask=[],extra_mask='',threshold=1.0,varianceindex=[],separate_variance_maps=[],doqu=False,minplots=False):
	print('Using Nside=' + str(nside))
	if len(rescale_amp) == 0:
		rescale_amp=np.ones(len(maps))
	if len(rescale_variance) == 0:
		rescale_variance=np.ones(len(maps))
	if len(apply_extra_mask) == 0:
		apply_extra_mask = np.zeros(len(maps))
	if extra_mask != '':
		extra_mask = hp.read_map(extra_mask)
	npix = hp.nside2npix(nside)

	if use_planck:
		planckmap = hp.read_map('/Users/mpeel/Documents/maps/wmap_planck_pol/weighted_both_debwk_feb2015_tqu.fits',field=None)
		# planckmap = hp.pixelfunc.reorder(planckmap[0],n2r=True)
		planckmap = hp.ud_grade(planckmap[0],nside,order_in='NEST',order_out='RING')

		planck_iqu = hp.read_map('/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits',field=None)
		# planckmap = hp.pixelfunc.reorder(planckmap[0],n2r=True)
		planck_iqu = hp.ud_grade(planck_iqu,nside,order_in='RING',order_out='RING')

	if use_cbass:
		# cbass = hp.read_map('/Users/mpeel/Documents/maps/cbass2019/512_60.00smoothed_cbass_4.76_512_mKCMBunits.fits',field=None)
		cbass = hp.read_map('/Users/mpeel/Documents/maps/cbass2019/cbass_global8p8deg_swapQU_NIGHT_v28allelsNs_37_noiseCut_masked5pc_G_1024_ol500_lessTol_g_map_g_1deg_0256.fits',field=None)
		# planckmap = hp.pixelfunc.reorder(planckmap[0],n2r=True)
		cbass = hp.ud_grade(cbass,nside,order_in='RING',order_out='RING')
		newmask = np.ones(npix)
		newmask[cbass[0] < -100.0] = 0
		hp.mollview(cbass[1]*1e3*((28.4/4.76)**index)*newmask,min=-0.05,max=0.05, title='CBASS Q')
		plt.savefig(outdirectory+'cbassq.png')
		plt.close()
		plt.clf()
		hp.mollview(cbass[2]*1e3*((28.4/4.76)**index)*newmask,min=-0.05,max=0.05,title='CBASS U')
		plt.savefig(outdirectory+'cbassu.png')
		plt.close()
		plt.clf()
		hp.mollview(planck_iqu[1]*newmask,min=-0.05,max=0.05,title='Planck Q')
		plt.savefig(outdirectory+'planckq.png')
		plt.close()
		plt.clf()
		hp.mollview(planck_iqu[2]*newmask,min=-0.05,max=0.05, title='Planck U')
		plt.savefig(outdirectory+'plancku.png')
		plt.close()
		plt.clf()
		print(max(planck_iqu[1]))
		print(max(cbass[1]))
		print(max(cbass[1])*(28.4/4.76)**index)
		cbass = cbass * 1e3
		if use_planck:
			# Output a comparison between C-BASS and the Planck data
			hp.mollview(((np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index)-1000.0*planckmap)*newmask,min=-0.05,max=0.05,title='CBASS - (Planck+WMAP)')
			plt.savefig(outdirectory+'cbass_diff_to_planckwmap.png')
			plt.close()
			plt.clf()
			hp.mollview((1000.0*planckmap-(np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index))*newmask,min=-0.05,max=0.05)
			plt.savefig(outdirectory+'cbass_diff_to_planckwmap_inverse.pdf')
			plt.close()
			plt.clf()

			hp.mollview(((np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index))*newmask,min=0,max=0.05,title='CBASS P')
			plt.savefig(outdirectory+'cbass_P.png')
			plt.close()
			plt.clf()
			hp.mollview((np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*newmask,min=0,max=0.05,title='Planck P')
			plt.savefig(outdirectory+'planck_P.png')
			plt.close()
			plt.clf()

			hp.mollview(((np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index)-np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*newmask,min=-0.03,max=0.03,title='CBASS - Planck')
			plt.savefig(outdirectory+'cbass_diff_to_planck.png')
			plt.close()
			plt.clf()

			hp.mollview((np.sqrt((cbass[1]*(28.4/4.76)**index - planck_iqu[1])**2+(cbass[2]*(28.4/4.76)**index - planck_iqu[2])**2))*newmask,max=0.05,title='CBASS - Planck via sqrt((QCB-QPlanck)**2 + (UCB-UPlanck)**2)')#,norm=colors.PowerNorm(gamma=0.2))
			plt.savefig(outdirectory+'cbass_diff_to_planckQU.png')
			plt.close()
			plt.clf()
			hp.mollview((cbass[1]*(28.4/4.76)**index - planck_iqu[1])*newmask,min=-0.05,max=0.05,title='CBASS Q - Planck Q')
			plt.savefig(outdirectory+'cbass_diff_to_planckQ.png')
			plt.close()
			plt.clf()

			hp.mollview((cbass[2]*(28.4/4.76)**index - planck_iqu[2])*newmask,min=-0.05,max=0.05,title='CBASS U - Planck U')
			plt.savefig(outdirectory+'cbass_diff_to_planckU.png')
			plt.close()
			plt.clf()

	nummaps = len(maps)
	commonmask = np.ones(npix)
	if doqu:
		combine = np.asarray([[np.zeros(npix), np.zeros(npix)],[np.zeros(npix), np.zeros(npix)]]).T
		weight = np.asarray([[np.zeros(npix), np.zeros(npix)],[np.zeros(npix), np.zeros(npix)]]).T
	else:
		combine_q = np.zeros(npix)
		combine_u = np.zeros(npix)
		weight_q = np.zeros(npix)
		weight_u = np.zeros(npix)
	# if doqu:
	# 	combine_qu = np.zeros(npix)
	# 	weight_qu = np.zeros(npix)

	rescale_vals = np.zeros((3,nummaps))
	for i in range(0,nummaps):
		print(maps[i])
		mapdata = hp.read_map(indirectory+maps[i],field=None)
		mapdata = hp.ud_grade(mapdata,nside,order_in='RING',order_out='RING')
		commonmask[mapdata[0][:] == hp.UNSEEN] = 0
		mapdata[0][mapdata[0][:] == hp.UNSEEN] = 0.0
		mapdata[1][mapdata[1][:] == hp.UNSEEN] = 0.0
		mapdata[2][mapdata[2][:] == hp.UNSEEN] = 0.0
		if not minplots:
			hp.mollview(mapdata[0],norm='hist')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_I.pdf')
			plt.close()
			plt.clf()
			hp.mollview(mapdata[1],norm='hist')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_Q.pdf')
			plt.close()
			plt.clf()
			hp.mollview(mapdata[2],norm='hist')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_U.pdf')
			plt.close()
			plt.clf()
			hp.mollview(np.sqrt(mapdata[1]**2+mapdata[2]**2),min=0,max=np.sqrt(threshold**2+threshold**2),cmap=plt.get_cmap('jet'))
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_P.pdf')
			plt.close()
			plt.clf()

		temp_Q = mapdata[1]*rescale_amp[i]*(normfreq/freqs[i])**index
		temp_U = mapdata[2]*rescale_amp[i]*(normfreq/freqs[i])**index
		if not minplots:
			hp.mollview(temp_Q,min=-threshold,max=threshold)
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_I_rescale.pdf')
			plt.close()
			plt.clf()
			hp.mollview(temp_U,min=-threshold,max=threshold)
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_Q_rescale.pdf')
			plt.close()
			plt.clf()
			hp.mollview(np.sqrt(temp_Q**2+temp_U**2),min=0,max=np.sqrt(threshold**2+threshold**2),cmap=plt.get_cmap('jet'))
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_U_rescale.pdf')
			plt.close()
			plt.clf()
		temp_Q = 0
		temp_U = 0


		if 'Planck' in maps[i]:

			if separate_variance_maps == [] and varianceindex == []:
				# var_i = hp.read_map(indirectory+maps[i].replace('tqu','tqu_noise').replace('512_','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_4_variance_'+str(nside)),field=None)
				var_q = hp.read_map(indirectory+maps[i].replace('tqu','tqu_noise').replace('512_','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_7_variance_'+str(nside)),field=None)
				var_u = hp.read_map(indirectory+maps[i].replace('tqu','tqu_noise').replace('512_','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_9_variance_'+str(nside)),field=None)

				# var_i = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_','').replace('256','1024').replace('bpcorr','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_4_actualvariance'),field=None)
				# var_q = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_','').replace('256','1024').replace('bpcorr','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_7_actualvariance'),field=None)
				# try:
				# 	var_u = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_','').replace('256','1024').replace('bpcorr','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_9_actualvariance'),field=None)
				# except:
				# 	var_u = var_q.copy()

				# var_i = hp.read_map(indirectory+maps[i].replace('2048_20.00s','20.0s').replace('_mKCMBunits','_mKCMBunits_noisenum2_4_actualvariance'),field=None)
				# var_q = hp.read_map(indirectory+maps[i].replace('2048_60.0s','60.0s').replace('20.00s','20.0s').replace('_mKCMBunits','_mKCMBunits_noisenum2_7_actualvariance'),field=None)
				# try:
				# 	var_u = hp.read_map(indirectory+maps[i].replace('2048_60.0s','60.0s').replace('20.00s','20.0s').replace('_mKCMBunits','_mKCMBunits_noisenum2_9_actualvariance'),field=None)
				# except:
				# 	var_u = var_q.copy()
			elif separate_variance_maps != []:
				# var_i = hp.read_map(indirectory+separate_variance_maps[i],field=None)[varianceindex[i][0]]
				var_q = hp.read_map(indirectory+separate_variance_maps[i],field=None)[varianceindex[i][1]]
				var_u = hp.read_map(indirectory+separate_variance_maps[i],field=None)[varianceindex[i][2]]
			else:
				# var_i = mapdata[varianceindex[i][0]].copy()
				var_q = mapdata[varianceindex[i][1]].copy()
				var_u = mapdata[varianceindex[i][2]].copy()
				if doqu:
					var_qu = mapdata[varianceindex[i][3]].copy()

			# var_i = hp.ud_grade(var_i,nside,power=2)
			var_q = hp.ud_grade(var_q,nside,power=2)
			var_u = hp.ud_grade(var_u,nside,power=2)
			if doqu:
				var_qu = hp.ud_grade(var_qu,nside,power=2)

			# var_i = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_','').replace('256','1024').replace('bpcorr','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_4_variance_'+str(nside)),field=None)
			# var_q = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_','').replace('256','1024').replace('bpcorr','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_7_variance_'+str(nside)),field=None)
			# try:
			# 	var_u = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_','').replace('256','1024').replace('bpcorr','').replace('60.0s','60.00s').replace('_mKCMBunits','_mKCMBunits_9_variance_'+str(nside)),field=None)
			# except:
			# 	var_u = var_q.copy()
			# var_i = mapdata[4].copy()
			# var_q = mapdata[7].copy()
			# var_u = mapdata[9].copy()
		elif 'wmap' in maps[i]:
			if varianceindex == []:
				# var_i = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_6','6').replace('60.0s','60.0s').replace('_mKCMBunits','_mKCMBunits_0.fits_actualvariance').replace('tqu','tqu_noise'),field=None)
				var_q = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_6','6').replace('60.0s','60.0s').replace('_mKCMBunits','_mKCMBunits_1.fits_actualvariance').replace('tqu','tqu_noise'),field=None)
				var_u = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_6','6').replace('60.0s','60.0s').replace('_mKCMBunits','_mKCMBunits_3.fits_actualvariance').replace('tqu','tqu_noise'),field=None)
				# var_i = hp.ud_grade(var_i,nside,power=2)
				var_q = hp.ud_grade(var_q,nside,power=2)
				var_u = hp.ud_grade(var_u,nside,power=2)
				# var_i = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_6','6').replace('60.0s','60.0s').replace('_mKCMBunits','_mKCMBunits_0.fits_variance_'+str(nside)).replace('tqu','tqu_noise'),field=None)
				# var_q = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_6','6').replace('60.0s','60.0s').replace('_mKCMBunits','_mKCMBunits_1.fits_variance_'+str(nside)).replace('tqu','tqu_noise'),field=None)
				# var_u = hp.read_map(indirectory.replace('tqu','tqu_noise')+maps[i].replace('512_6','6').replace('60.0s','60.0s').replace('_mKCMBunits','_mKCMBunits_3.fits_variance_'+str(nside)).replace('tqu','tqu_noise'),field=None)
				# var_i = mapdata[3].copy()
				# var_q = mapdata[3].copy()
				# var_u = mapdata[3].copy()
			else:
				# var_i = mapdata[varianceindex[i][0]].copy()
				var_q = mapdata[varianceindex[i][1]].copy()
				var_u = mapdata[varianceindex[i][2]].copy()
				if doqu:
					var_qu = mapdata[varianceindex[i][3]].copy()

			# var_i = hp.ud_grade(var_i,nside,power=2)
			var_q = hp.ud_grade(var_q,nside,power=2)
			var_u = hp.ud_grade(var_u,nside,power=2)
			if doqu:
				var_qu = hp.ud_grade(var_qu,nside,power=2)

		if use_halfrings:
			map_half1 = hp.read_map(indirectory+maps_half1[i],field=None)
			map_half2 = hp.read_map(indirectory+maps_half2[i],field=None)
			# var_i = (np.abs(map_half1[0] - map_half2[0])/2.0)**2
			var_q = (np.abs(map_half1[1] - map_half2[1])/2.0)**2
			var_u = (np.abs(map_half1[2] - map_half2[2])/2.0)**2
			# var_i[var_i == 0.0] = 1e4
			var_q[var_q == 0.0] = 1e4
			var_u[var_u == 0.0] = 1e4
			# var_i[var_i < np.median(var_i)] = np.median(var_i)
			var_q[var_q < np.median(var_q)] = np.median(var_q)
			var_u[var_u < np.median(var_u)] = np.median(var_u)
			# var_i[var_i > 1e4] = 1e4
			var_q[var_q > 1e4] = 1e4
			var_u[var_u > 1e4] = 1e4
		elif use_weights:
			# Get the variance maps
			# var_i = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_0_variance'),field=None)
			var_q = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
			var_u = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
		elif use_reweight_by_rms and 'Planck' not in maps[i] and 'wmap' not in maps[i]:
			map_half1 = hp.read_map(indirectory+maps_half1[i],field=None)
			map_half2 = hp.read_map(indirectory+maps_half2[i],field=None)
			# diff_i = (np.abs(map_half1[0] - map_half2[0])/2.0)#**2
			diff_q = (np.abs(map_half1[1] - map_half2[1])/2.0)#**2
			diff_u = (np.abs(map_half1[2] - map_half2[2])/2.0)#**2
			# diff_i[diff_i > 1e4] = 0.0
			diff_q[diff_q > 1e4] = 0.0
			diff_u[diff_u > 1e4] = 0.0
			# print(np.std(diff_i[diff_i != 0.0]))
			# print(np.std(diff_q[diff_q != 0.0]))
			# print(np.std(diff_u[diff_u != 0.0]))
			# var_i = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_0_variance'),field=None)
			var_q = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
			var_u = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
			# var_i[var_i < -1e4] = 0.0
			var_q[var_q < -1e4] = 0.0
			var_u[var_u < -1e4] = 0.0
			# print(np.max(var_i))
			# print(np.min(var_i))

			if use_reweight_by_rms_method == 1:
				# noise_i = noiserealisation(np.sqrt(var_i[var_i != 0.0]),len(var_i[var_i != 0.0]))
				noise_q = noiserealisation(np.sqrt(var_q[var_q != 0.0]),len(var_q[var_q != 0.0]))
				noise_u = noiserealisation(np.sqrt(var_u[var_u != 0.0]),len(var_u[var_u != 0.0]))
				rescale_vals[0,i] = np.std(diff_i[diff_i != 0.0])/np.std(noise_i[noise_i != 0.0])
				rescale_vals[1,i] = np.std(diff_q[diff_q != 0.0])/np.std(noise_q[noise_q != 0.0])
				rescale_vals[2,i] = np.std(diff_u[diff_u != 0.0])/np.std(noise_u[noise_u != 0.0])
			else:
				# rescale_vals[0,i] = np.std(diff_i[var_i != 0.0] / np.sqrt(var_i[var_i != 0.0]))
				rescale_vals[1,i] = np.std(diff_q[var_q != 0.0] / np.sqrt(var_q[var_q != 0.0]))
				rescale_vals[2,i] = np.std(diff_u[var_u != 0.0] / np.sqrt(var_u[var_u != 0.0]))

			# var_i[:] = var_i[:] * (rescale_vals[0,i])**2.0
			var_q[:] = var_q[:] * (rescale_vals[1,i])**2.0
			var_u[:] = var_u[:] * (rescale_vals[2,i])**2.0
			# var_i[var_i == 0.0] = 1e4
			var_q[var_q == 0.0] = 1e4
			var_u[var_u == 0.0] = 1e4


		# print(maps[i])
		# print(np.median(np.sqrt(var_i[var_i[:] >=0])))
		# print(np.median(np.sqrt(var_q[var_q[:] >=0])))
		# print(np.median(np.sqrt(var_u[var_u[:] >=0])))

		if not minplots:
			hp.mollview(var_q,norm='hist')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_Q_var.pdf')
			plt.close()
			plt.clf()
			hp.mollview(var_u,norm='hist')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_U_var.pdf')
			plt.close()
			plt.clf()
			if doqu:
				hp.mollview(var_qu,norm='hist')
				plt.savefig(outdirectory+maps[i].split('/')[-1]+'_QU_var.pdf')
				plt.close()
				plt.clf()

		# var_i = var_i * ((normfreq/freqs[i])**index)**2
		var_q = var_q * ((normfreq/freqs[i])**index)**2
		var_u = var_u * ((normfreq/freqs[i])**index)**2
		if doqu:
			var_qu = var_qu * ((normfreq/freqs[i])**index)**2


		# hp.mollview(var_i,norm='hist')
		# plt.savefig(outdirectory+maps[i].split('/')[-1]+'_0_var.pdf')
		if i == 0:
			qmin = np.min(var_q*rescale_variance[i])
			qmax = np.max(var_q*rescale_variance[i])/4.0
			umin = np.min(var_u*rescale_variance[i])
			umax = np.max(var_u*rescale_variance[i])/4.0
		if not minplots:
			hp.mollview(var_q*rescale_variance[i],min=qmin*(normfreq/freqs[i])**index,max=qmax*(normfreq/freqs[i])**index)#,norm='hist')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_Q_var_rescale.pdf')
			plt.close()
			plt.clf()
			hp.mollview(var_u*rescale_variance[i],min=umin*(normfreq/freqs[i])**index,max=umax*(normfreq/freqs[i])**index)#,norm='hist')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'_U_var_rescale.pdf')
			plt.close()
			plt.clf()
			if doqu:
				hp.mollview(var_qu*rescale_variance[i],min=-umin*(normfreq/freqs[i])**index,max=umax*(normfreq/freqs[i])**index)#,norm='hist')
				plt.savefig(outdirectory+maps[i].split('/')[-1]+'_QU_var_rescale.pdf')
				plt.close()
				plt.clf()

		# print(np.median(np.sqrt(var_i[var_i[:] >=0])))
		# print(np.median(np.sqrt(var_q[var_q[:] >=0])))
		# print(np.median(np.sqrt(var_u[var_u[:] >=0])))

		# Output a comparison with the Planck data
		if use_planck:
			hp.mollview(((np.sqrt(mapdata[1]**2+mapdata[2]**2)*(28.4/freqs[i])**index)-1000.0*planckmap)*commonmask,min=-0.1,max=0.1,title=maps[i].split('/')[-1] + ' - (Planck+WMAP)')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckwmap.png')
			plt.close()
			plt.clf()
			hp.mollview((1000.0*planckmap-(np.sqrt(mapdata[1]**2+mapdata[2]**2)*(28.4/freqs[i])**index))*commonmask,min=-0.1,max=0.1)
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckwmap_inverse.pdf')
			plt.close()
			plt.clf()

			hp.mollview((np.sqrt((mapdata[1]*(28.4/freqs[i])**index)**2+(mapdata[2]*(28.4/freqs[i])**index)**2)-np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*commonmask,min=-0.05,max=0.05,title=maps[i].split('/')[-1] + ' P - Planck P')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckP.png')
			plt.close()
			plt.clf()
			hp.mollview(np.sqrt((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])**2+(mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])**2)*commonmask,min=0,max=0.1,title=maps[i].split('/')[-1] + ' - Planck via sqrt((QMFI-QPlanck)**2 + (UMFI-UPlanck)**2)')#,norm=colors.PowerNorm(gamma=0.2))
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckQU.png')
			plt.close()
			plt.clf()
			hp.mollview((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])*commonmask,min=-0.05,max=0.05,title=maps[i].split('/')[-1] + ' Q - Planck Q')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckQ.png')
			plt.close()
			plt.clf()

			hp.mollview((mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])*commonmask,min=-0.05,max=0.05,title=maps[i].split('/')[-1] + ' U - Planck U')
			plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckU.png')
			plt.close()
			plt.clf()

			# hp.gnomview(((mapdata[0]*(28.4/freqs[i])**index))*commonmask,title=maps[i].split('/')[-1] + ' I (rescaled to 28.4)',rot=[80,0],max=10.0,reso=10.0)
			# plt.savefig(outdirectory+maps[i].split('/')[-1]+'_I_cyg.png')
			# plt.close()
			# plt.clf()
			# hp.gnomview(((mapdata[0]*(28.4/freqs[i])**index)-planck_iqu[0])*commonmask,min=-2.0,max=2.0,title=maps[i].split('/')[-1] + ' I - Planck I',rot=[80,0],reso=10.0)
			# plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckI_cyg.png')
			# plt.close()
			# plt.clf()
			# hp.gnomview((np.sqrt((mapdata[1]*(28.4/freqs[i])**index)**2+(mapdata[2]*(28.4/freqs[i])**index)**2)-np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*commonmask,min=-0.05,max=0.05,rot=[80,0],reso=10.0,title=maps[i].split('/')[-1] + ' P - Planck P')
			# plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckP_cyg.png')
			# plt.close()
			# plt.clf()
			# hp.gnomview(np.sqrt((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])**2+(mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])**2)*commonmask,min=0,max=0.1,rot=[80,0],reso=10.0,title=maps[i].split('/')[-1] + ' - Planck via sqrt((QMFI-QPlanck)**2 + (UMFI-UPlanck)**2)')#,norm=colors.PowerNorm(gamma=0.2))
			# plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckQU_cyg.png')
			# plt.close()
			# plt.clf()
			# hp.gnomview((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])*commonmask,min=-0.1,max=0.1,rot=[80,0],reso=10.0,title=maps[i].split('/')[-1] + ' Q - Planck Q')
			# plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckQ_cyg.png')
			# plt.close()
			# plt.clf()

			# hp.gnomview((mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])*commonmask,min=-0.1,max=0.1,rot=[80,0],reso=10.0,title=maps[i].split('/')[-1] + ' U - Planck U')
			# plt.savefig(outdirectory+maps[i].split('/')[-1]+'diff_to_planckU_cyg.png')
			# plt.close()
			# plt.clf()
			# exit()


		if apply_extra_mask[i]:
			# Scale up variances to get rid of those values
			var_q[extra_mask == 0] *= 1e30
			var_u[extra_mask == 0] *= 1e30
			if doqu:
				var_qu[extra_mask == 0] *= 1e30


		snmap = (np.sqrt(mapdata[1]**2+mapdata[2]**2)*rescale_amp[i]*(normfreq/freqs[i])**index)*commonmask/np.sqrt((var_q*rescale_variance[i])+(var_u*rescale_variance[i]))

		print('Median S/N:' + str(np.median(snmap[commonmask==1])))
		print('Mean S/N:' + str(np.mean(snmap[commonmask==1])))

		p_temp = np.sqrt((mapdata[1]*rescale_amp[i]*(normfreq/freqs[i])**index)**2+(mapdata[2]*rescale_amp[i]*(normfreq/freqs[i])**index)**2)
		var_p = ((var_q*(mapdata[1]*rescale_amp[i]*(normfreq/freqs[i])**index)**2)+(var_u*(mapdata[2]*rescale_amp[i]*(normfreq/freqs[i])**index)**2))/(p_temp*p_temp)

		snmap = p_temp/np.sqrt(var_p)

		print('Median S/N:' + str(np.median(snmap[commonmask==1])))
		print('Mean S/N:' + str(np.mean(snmap[commonmask==1])))
		# exit()

		if doqu:
			W = np.linalg.inv(np.asarray([[var_q*rescale_variance[i],var_qu*rescale_variance[i]],[var_qu*rescale_variance[i],var_u*rescale_variance[i]]]).T)
			D = np.asarray([[mapdata[1]*rescale_amp[i]*(normfreq/freqs[i])**index,np.zeros(len(mapdata[1]))],[np.zeros(len(mapdata[1])),mapdata[2]*rescale_amp[i]*(normfreq/freqs[i])**index]]).T
			combine = combine + D @ W
			weight = weight + W
		else:
			if i == 0:
				combine_q = (mapdata[1].copy()*rescale_amp[i]*(normfreq/freqs[i])**index)/(var_q*rescale_variance[i])
				combine_u = (mapdata[2].copy()*rescale_amp[i]*(normfreq/freqs[i])**index)/(var_u*rescale_variance[i])
				weight_q = 1.0/(var_q.copy()*rescale_variance[i])
				weight_u = 1.0/(var_u.copy()*rescale_variance[i])
			else:
				combine_q = combine_q+(mapdata[1]*rescale_amp[i]*(normfreq/freqs[i])**index)/(var_q*rescale_variance[i])
				combine_u = combine_u+(mapdata[2]*rescale_amp[i]*(normfreq/freqs[i])**index)/(var_u*rescale_variance[i])
				weight_q = weight_q+1.0/(var_q*rescale_variance[i])
				weight_u = weight_u+1.0/(var_u*rescale_variance[i])
	

	if doqu:
		print(combine)
		print(np.max(combine))
		print(np.shape(combine))
		print(weight)
		print(np.shape(weight))
		# for l in range(0,len(combine[0])):
		# 	combine[l] = combine[l] @ np.linalg.inv(weight[l])
		inv_W = np.linalg.inv(weight)
		combine = combine @ inv_W
		print(combine)
		print(np.shape(combine))
		print(np.max(combine))
		# exit()
		combine_q = combine[:,0,0]
		combine_u = combine[:,1,1]
		combine_qu = combine[:,0,1]
		weight_q = 1.0/inv_W[:,0,0]
		weight_u = 1.0/inv_W[:,1,1]
		weight_qu = 1.0/inv_W[:,0,1]
	else:
		combine_q /= weight_q
		combine_u /= weight_u
	
	if not minplots:
		hp.write_map(outdirectory+prefix+'_commonmask.fits',commonmask,overwrite=True)
		hp.mollview(commonmask)
		plt.savefig(outdirectory+prefix+'_commonmask.pdf')
		plt.close()
		plt.clf()

		hp.write_map(outdirectory+prefix+'_combine_q.fits',combine_q*commonmask,overwrite=True)
		hp.write_map(outdirectory+prefix+'_combine_q_unc.fits',(1.0/weight_q)*commonmask,overwrite=True)
		hp.mollview(combine_q*commonmask,min=-threshold,max=threshold)
		plt.savefig(outdirectory+prefix+'_combine_q.pdf')
		plt.close()
		plt.clf()
		hp.mollview(combine_qu*commonmask,min=-threshold,max=threshold)
		plt.savefig(outdirectory+prefix+'_combine_qu.pdf')
		plt.close()
		plt.clf()

		hp.write_map(outdirectory+prefix+'_combine_u.fits',combine_u*commonmask,overwrite=True)
		hp.write_map(outdirectory+prefix+'_combine_u_unc.fits',(1.0/weight_u)*commonmask,overwrite=True)
		hp.mollview(combine_u*commonmask,min=-threshold,max=threshold)
		plt.savefig(outdirectory+prefix+'_combine_u.pdf')
		plt.close()
		plt.clf()

		hp.write_map(outdirectory+prefix+'_combine_P.fits',np.sqrt(combine_q**2+combine_u**2)*commonmask,overwrite=True)
		hp.mollview(np.sqrt(combine_q**2+combine_u**2)*commonmask,min=0,max=np.sqrt(threshold**2+threshold**2),cmap=plt.get_cmap('jet'))
		plt.savefig(outdirectory+prefix+'_combine_P.pdf')
		plt.close()
		plt.clf()

	hp.write_map(outdirectory+prefix+'_combine_P_nomask.fits',np.sqrt(combine_q**2+combine_u**2),overwrite=True)
	hp.mollview(np.sqrt(combine_q**2+combine_u**2),min=0,max=np.sqrt(threshold**2+threshold**2),cmap=plt.get_cmap('jet'))
	plt.savefig(outdirectory+prefix+'_combine_P_nomask.pdf')
	plt.close()
	plt.clf()


	cols = []
	cols.append(fits.Column(name='P', format='E', array=np.asarray(np.sqrt(combine_q**2+combine_u**2))))
	cols.append(fits.Column(name='Q', format='E', array=np.asarray(combine_q)))
	if doqu:
		cols.append(fits.Column(name='QU', format='E', array=np.asarray(combine_qu)))
	cols.append(fits.Column(name='U', format='E', array=np.asarray(combine_u)))
	cols.append(fits.Column(name='QQ_cov', format='E', array=np.asarray(1.0/weight_q)))
	if doqu:
		cols.append(fits.Column(name='QU_cov', format='E', array=np.asarray(1.0/weight_qu)))
	cols.append(fits.Column(name='UU_cov', format='E', array=np.asarray(1.0/weight_u)))

	cols = fits.ColDefs(cols)
	bin_hdu = fits.BinTableHDU.from_columns(cols)
	bin_hdu.header['ORDERING']='RING'
	bin_hdu.header['POLCONV']= 'COSMO'
	bin_hdu.header['PIXTYPE']= 'HEALPIX'
	bin_hdu.header['NSIDE']  = str(nside)
	bin_hdu.header['COMMENT']= "Weighted map calculated by Mike Peel's weighted_pol_map"
	bin_hdu.header['FREQ']= str(normfreq)
	bin_hdu.header['INDEX']= str(index)
	bin_hdu.header['INPMAPS']= str(maps)
	bin_hdu.header['TUNIT1'] = 'mK_CMB'
	bin_hdu.header['TUNIT2'] = 'mK_CMB'
	bin_hdu.header['TUNIT3'] = 'mK_CMB'
	if doqu:
		bin_hdu.header['TUNIT4'] = 'mK_CMB'
		bin_hdu.header['TUNIT5'] = '(mK_CMB)^2'
		bin_hdu.header['TUNIT6'] = '(mK_CMB)^2'
		bin_hdu.header['TUNIT7'] = '(mK_CMB)^2'
		bin_hdu.header['TUNIT8'] = '(mK_CMB)^2'
	else:
		bin_hdu.header['TUNIT4'] = '(mK_CMB)^2'
		bin_hdu.header['TUNIT5'] = '(mK_CMB)^2'
		bin_hdu.header['TUNIT6'] = '(mK_CMB)^2'
	# Write out the file
	bin_hdu.writeto(outdirectory+prefix+'_combine.fits',overwrite=True)

	# hp.write_map(outdirectory+prefix+'_combine.fits',[np.sqrt(combine_q**2+combine_u**2),combine_q,combine_u,1.0/weight_q,1.0/weight_u],overwrite=True)

	if not minplots:
		hp.mollview(1.0/weight_q,min=qmin,max=qmax)#,norm='hist')
		plt.savefig(outdirectory+prefix+'_combine_Q_var.pdf')
		plt.close()
		plt.clf()
		hp.mollview(1.0/weight_u,min=umin,max=umax)#,norm='hist')
		plt.savefig(outdirectory+prefix+'_combine_U_var.pdf')
		plt.close()
		plt.clf()
		if doqu:
			hp.mollview(1.0/weight_qu,min=umin,max=umax)#,norm='hist')
			plt.savefig(outdirectory+prefix+'_combine_QU_var.pdf')
			plt.close()
			plt.clf()

		snmap = np.sqrt(combine_q**2+combine_u**2)*commonmask/np.sqrt(1.0/weight_q+1.0/weight_u)

		print('Median S/N:' + str(np.median(snmap[commonmask==1])))
		print('Mean S/N:' + str(np.mean(snmap[commonmask==1])))
		hp.mollview(snmap,min=0,max=3.0)#,norm='hist')
		plt.savefig(outdirectory+prefix+'_snmap.pdf')
		plt.close()
		plt.clf()

	# commonmask2 = hp.ud_grade(commonmask,256,order_in='RING',order_out='RING')
	if use_planck:
		hp.write_map(outdirectory+'wmapplanck2015.fits',[np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2),planck_iqu[1],planck_iqu[2]],overwrite=True)
		hp.mollview(planckmap*1000.0*commonmask,min=0,max=np.sqrt(threshold**2+threshold**2),cmap=plt.get_cmap('jet'))
		plt.savefig(outdirectory+'combine_P_planck.pdf')
		plt.close()
		plt.clf()
		hp.mollview(planckmap*1000.0,min=0,max=0.06,cmap=plt.get_cmap('jet'))
		plt.savefig(outdirectory+'combine_P_planck_nomask.pdf')
		plt.close()
		plt.clf()

	if use_reweight_by_rms:
		np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
		print(rescale_vals)

	return