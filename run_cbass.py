import astropy.io.fits as fits
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import os

from smoothmap import *

output_resolution = 60.0
output_nside = [512]#, 256, 128, 64, 32, 16, 8]
smoothvariance = False

directory = '/Users/mpeel/Documents/maps/'
outdirectory = directory+"cbass2020_tqu_v1.4/"
os.makedirs(outdirectory, exist_ok=True)

beamtf_p30 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',28)
print(beamtf_p30)

beam_cbass_angle_first, beam_cbass_val_first, beam_int_first = np.loadtxt(directory+'cbass2020/cbass-beam.txt',unpack=True)
beam_cbass_val_first /= beam_cbass_val_first[0]
plt.plot(beam_cbass_angle_first*180.0/np.pi,beam_cbass_val_first)
beam_cbass_angle, beam_cbass_val = np.loadtxt(directory+'cbass2020/cbass_beam2.txt',unpack=True)
print(beam_cbass_val)
plt.plot(beam_cbass_angle,beam_cbass_val)
plt.yscale('log')
plt.savefig(outdirectory+'beams_orig.pdf')
plt.clf()

cbass_wf = hp.sphtfunc.beam2bl(beam_cbass_val,beam_cbass_angle*np.pi/180.0,3*output_nside[0])
cbass_wf -= np.min(cbass_wf)
cbass_wf /= cbass_wf[0]
print(np.min(cbass_wf))
cbass_wf_first = hp.sphtfunc.beam2bl(beam_cbass_val_first,beam_cbass_angle_first,3*output_nside[0])
cbass_wf_first -= np.min(cbass_wf_first)
cbass_wf_first /= cbass_wf_first[0]
print(np.min(cbass_wf))
# cbass_wf = cbass_wf**0.5
print(cbass_wf)
plt.plot(cbass_wf,label='cbass')
plt.plot(cbass_wf_first,label='cbass_first')
plt.plot(beamtf_p30,label='planck')
plt.yscale('log')

cbass_tf_I_l, cbass_tf_I_val, cbass_tf_I_gauss, cbass_tf_I_diff = np.loadtxt(directory+'cbass2020/smoothKink_lmin555_lmax800_windowFunctionsI.txt')
plt.plot(cbass_tf_I_val,label='cbass_tf')
plt.legend()
plt.xlim(0,1000)

plt.savefig(outdirectory+'beams.pdf')
plt.clf()
# exit()
numnside = len(output_nside)
for i in range(0,numnside):
	# smoothmap(directory,outdirectory,'cbass2020/v33b/NIGHTMERID20/calibrated_map/tauA_cal_NM20_v33b_allelsNS1_xAS14_masked5pc_G_1024_ol500_lessTol_g_Pipe_map.fits',str(output_nside[i])+'_60.00smoothed_cbassv33bBeam_4.76_512_mKCMBunits.fits', output_resolution, nside_out=output_nside[i], windowfunction=cbass_tf_I_val,usehealpixfits=True,minmapvalue=-1e10,maxmapvalue=1e10,minmaxmaps=[0,1,2],smoothvariance=smoothvariance)
	# smoothmap(directory,outdirectory,'cbass2020/v33b/NIGHTMERID20/calibrated_map/tauA_cal_NM20_v33b_allelsNS1_xAS14_masked5pc_G_1024_ol500_lessTol_g_Pipe_map.fits',str(output_nside[i])+'_60.00smoothed_cbassv33bBeam2_4.76_512_mKCMBunits.fits', output_resolution, nside_out=output_nside[i], windowfunction=cbass_wf,usehealpixfits=True,minmapvalue=-1e10,maxmapvalue=1e10,minmaxmaps=[0,1,2],smoothvariance=smoothvariance,taper=True,lmin_taper=450,lmax_taper=600)
	# smoothmap(directory,outdirectory,'cbass2020/v33b/NIGHTMERID20/calibrated_map/tauA_cal_NM20_v33b_allelsNS1_xAS14_masked5pc_G_1024_ol500_lessTol_g_Pipe_map.fits',str(output_nside[i])+'_60.00smoothed_cbassv33bBeam3_4.76_512_mKCMBunits.fits', output_resolution, nside_out=output_nside[i], windowfunction=cbass_wf,usehealpixfits=True,minmapvalue=-1e10,maxmapvalue=1e10,minmaxmaps=[0,1,2],smoothvariance=smoothvariance,taper_gauss=True,taper_gauss_sigma=45.0)
	# smoothmap(directory,outdirectory,'cbass2020/v33b/NIGHTMERID20/calibrated_map/tauA_cal_NM20_v33b_allelsNS1_xAS14_masked5pc_G_1024_ol500_lessTol_g_Pipe_map.fits',str(output_nside[i])+'_45.00smoothed_cbassv33bBeam3_4.76_512_mKCMBunits.fits', 45.0, nside_out=output_nside[i], windowfunction=cbass_wf,usehealpixfits=True,minmapvalue=-1e10,maxmapvalue=1e10,minmaxmaps=[0,1,2],smoothvariance=smoothvariance,taper_gauss=True,taper_gauss_sigma=45.0)
	smoothmap(directory,outdirectory,'cbass2020/v33b/NIGHTMERID20/calibrated_map/tauA_cal_NM20_v33b_allelsNS1_xAS14_masked5pc_G_1024_ol500_lessTol_g_Pipe_map.fits',str(output_nside[i])+'_60.00smoothed_cbassv33bBeam4_4.76_512_mKCMBunits.fits', 60.0, nside_out=output_nside[i], windowfunction=cbass_wf_first,usehealpixfits=True,minmapvalue=-1e10,maxmapvalue=1e10,minmaxmaps=[0,1,2],smoothvariance=smoothvariance,taper=True,lmin_taper=450,lmax_taper=600)

	# smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_100_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_100_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p100,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)

# EOF
