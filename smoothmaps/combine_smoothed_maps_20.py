# Read in the various computed smoothed amplitude and noise maps, and save them into easier to use files
# Mike Peel		21 Sep 2020		Started
import numpy as np
import healpy as hp
import os
import astropy.io.fits as fits

def get_header_val(hdr,search):
	for i in range(0,len(hdr)):
		if search in hdr[i][0]:
			return hdr[i][1]
	return ''

def docombine(outfile, iqu_file, II_file, QQ_file='', UU_file='',rescale=1.0,comment='', QQ_add='', UU_add='',add_rescale=1.0):
	print(outfile)
	if os.path.isfile(outfile):
		print("You already have a file with the output name " + outfile + "! Not going to overwrite it. Move it, or set a new output filename, and try again!")
		return
	try:
		iqu,h = hp.read_map(iqu_file,field=None,h=True)
	except:
		iqu,h = hp.read_map(iqu_file.replace('20.0','20.00'),field=None,h=True)
	try:
		ii,hii = hp.read_map(II_file,h=True)
	except:
		ii,hii = hp.read_map(II_file.replace('20.0','20.00'),h=True)
	if QQ_file != '':
		try:
			qq,hqq = hp.read_map(QQ_file,field=None,h=True)
		except:
			qq,hqq = hp.read_map(QQ_file.replace('20.0','20.00'),h=True)
		try:
			uu,huu = hp.read_map(UU_file,h=True)
		except:
			uu,huu = hp.read_map(UU_file.replace('20.0','20.00'),h=True)
		try:
			qu,hqu = hp.read_map(QU_file,h=True)
		except:
			qu,hqu = hp.read_map(QU_file.replace('20.0','20.00'),h=True)
		print(np.shape(iqu[0]))
		print(np.shape(iqu[1]))
		print(np.shape(iqu[2]))
		print(np.shape(ii))
		print(np.shape(qq))
		print(np.shape(uu))
		print(np.shape(qu))
		if QQ_add != '':
			try:
				qq_add,hqq_add = hp.read_map(QQ_add,field=None,h=True)
			except:
				qq_add,hqq_add = hp.read_map(QQ_add.replace('20.0','20.00'),h=True)
			try:
				uu_add,huu_add = hp.read_map(UU_add,h=True)
			except:
				uu_add,huu_add = hp.read_map(UU_add.replace('20.0','20.00'),h=True)
			qq = qq + qq_add*add_rescale
			uu = uu + qq_add*add_rescale

	cols = []
	if QQ_file != '':
		cols.append(fits.Column(name='I', format='E', array=np.asarray(iqu[0])))
		cols.append(fits.Column(name='Q', format='E', array=np.asarray(iqu[1])))
		cols.append(fits.Column(name='U', format='E', array=np.asarray(iqu[2])))
		cols.append(fits.Column(name='II_cov', format='E', array=np.asarray(ii)*rescale))
		cols.append(fits.Column(name='QQ_cov', format='E', array=np.asarray(qq)*rescale))
		cols.append(fits.Column(name='QU_cov', format='E', array=np.asarray(qu)*rescale))
		cols.append(fits.Column(name='UU_cov', format='E', array=np.asarray(uu)*rescale))
	else:
		cols.append(fits.Column(name='I', format='E', array=np.asarray(iqu[0])))
		cols.append(fits.Column(name='II_cov', format='E', array=np.asarray(ii)*rescale))

	cols = fits.ColDefs(cols)
	bin_hdu = fits.BinTableHDU.from_columns(cols)
	bin_hdu.header['ORDERING']='RING'
	bin_hdu.header['POLCONV']='COSMO'
	bin_hdu.header['PIXTYPE']='HEALPIX'
	bin_hdu.header['COMMENT']=comment
	bin_hdu.header['NSIDE'] = get_header_val(hii,'NSIDE')
	bin_hdu.header['TUNIT1'] = get_header_val(h,'TUNIT1')
	if QQ_file != '':
		bin_hdu.header['TUNIT2'] = get_header_val(h,'TUNIT2')
		bin_hdu.header['TUNIT3'] = get_header_val(h,'TUNIT3')
		bin_hdu.header['TUNIT4'] = get_header_val(hii,'TUNIT1')
		bin_hdu.header['TUNIT5'] = get_header_val(hqq,'TUNIT1')
		bin_hdu.header['TUNIT6'] = get_header_val(hqu,'TUNIT1')
		bin_hdu.header['TUNIT7'] = get_header_val(huu,'TUNIT1')
	else:
		bin_hdu.header['TUNIT2'] = get_header_val(hii,'TUNIT1')
	bin_hdu.writeto(outfile)
	return 0

# directory = '/scratch1/mpeel/maps/'
directory = '/share/nas_cbassarc/mpeel/'
output_nside = [2048, 1024, 512, 256, 128, 64, 32, 16, 8]
comment = "Smoothed using Mike Peel's smoothmap.py v1.5 and smoothnoisemap.py v1.0"

# Planck 2020
mapdir = directory+'planck2020_tqu_v1.5/'
noisedir = directory+'planck2020_tqu_noise_v1.0/'
outdirectory = directory+"planck2020_tqu_v1.5_noise_v0.9_20arcmin/"
os.makedirs(outdirectory, exist_ok=True)
rescale=1000.0**2
for nside in output_nside:
	# Standard maps
	namestrings = ['100_2048_2020','143_2048_2020','217_2048_2020','353_2048_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_20.0smoothed_PlanckR4fullbeamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Dipole subtracted maps - using the same noise realisations
	namestrings = ['100_2048_2020','143_2048_2020','217_2048_2020','353_2048_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_20.0smoothed_PlanckR4fullbeamnodpNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_20.0smoothed_PlanckR4fullbeamnodp_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Intensity only maps
	namestrings = ['545_2048_2020','857_2048_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_20.0smoothed_PlanckR4fullbeamNoise_'+namestrings[i]+'_MJySrunits.fits',\
				mapdir+str(nside)+'_20.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_MJySrunits.fits',\
				noisedir+'20.0smoothed_PlanckR3fullbeam_'+namestrings[i]+'_MJySrunits_variance_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue

