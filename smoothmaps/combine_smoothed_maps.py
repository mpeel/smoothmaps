# Read in the various computed smoothed amplitude and noise maps, and save them into easier to use files
# Mike Peel		21 Sep 2020		Started
import numpy as np
import healpy as hp
import os

def docombine(outfile, iqu_file, II_file, QQ_file, UU_file,comment=''):
	iqu,h = hp.read_map(iqu_file,field=None,h=True)
	ii = hp.read_map(II_file,field=None)
	qq = hp.read_map(QQ_file,field=None)
	uu = hp.read_map(UU_file,field=None)

	cols = []
	cols.append(fits.Column(name='I', format='E', array=iqu[0]))
	cols.append(fits.Column(name='Q', format='E', array=iqu[1]))
	cols.append(fits.Column(name='U', format='E', array=iqu[2]))
	cols.append(fits.Column(name='II_cov', format='E', array=ii))
	cols.append(fits.Column(name='QQ_cov', format='E', array=qq))
	cols.append(fits.Column(name='UU_cov', format='E', array=uu))
	cols = fits.ColDefs(cols)
	bin_hdu = fits.BinTableHDU.from_columns(cols)
	bin_hdu.header['ORDERING']='RING'
	bin_hdu.header['POLCONV']='COSMO'
	bin_hdu.header['PIXTYPE']='HEALPIX'
	bin_hdu.header['COMMENT']=comment
	bin_hdu.writeto(outfile)
	return 0

# directory = '/scratch1/mpeel/maps/'
directory = '/share/nas_cbassarc/mpeel/'
output_nside = [2048, 1024, 512, 256, 128, 64, 32, 16, 8]
comment = "Smoothed using Mike Peel's smoothmap.py v1.4 and smoothnoisemap.py v0.7"
# WMAP9
mapdir = directory+'wmap9_tqu_v1.4'
noisedir = directory+'wmap9_tqu_noise_v0.7'
outdirectory = directory+"wmap9_tqu_v1.4_noise_v0.7/"
os.makedirs(outdirectory, exist_ok=True)

for nside in output_nside:
	if nside <= 512:
		namestrings = ['22.8_512_2013','33.0_512_2013','40.7_512_2013','60.7_512_2013','93.5_512_2013']
		for namestr in namestrings:
			docombine(outdirectory+str(nside)+'_60.0smoothed_wmap9beamNoise_'+namestr+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_wmap9beam_'+namestr+'_mKCMBunits.fits',[0,1,2],\
				noisedir+'60.00smoothed_wmap9beam_'+namestr+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.00smoothed_wmap9beam_'+namestr+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.00smoothed_wmap9beam_'+namestr+'_mKCMBunits_variance_U_'+str(nside)+'.fits',comment=comment)

# Planck 2018
# outdirectory = directory+"planck2018_tqu_v1.4_noise_v0.6/"
# os.makedirs(outdirectory, exist_ok=True)

# for nside in output_nside:
