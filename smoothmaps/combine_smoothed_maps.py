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

def docombine(outfile, iqu_file, II_file, QQ_file='', UU_file='',QU_file='',rescale=1.0,comment='', QQ_add='', UU_add='', QU_add='',add_rescale=1.0):
	print(outfile)
	if os.path.isfile(outfile):
		print("You already have a file with the output name " + outfile + "! Not going to overwrite it. Move it, or set a new output filename, and try again!")
		return
	try:
		iqu,h = hp.read_map(iqu_file,field=None,h=True)
		print(iqu_file)
	except:
		iqu,h = hp.read_map(iqu_file.replace('60.0','60.00'),field=None,h=True)
		print('2')
	try:
		ii,hii = hp.read_map(II_file,h=True)
		print(II_file)
	except:
		ii,hii = hp.read_map(II_file.replace('60.0','60.00'),h=True)
		print('3')
	if QQ_file != '':
		try:
			qq,hqq = hp.read_map(QQ_file,field=None,h=True)
			print(QQ_file)
		except:
			qq,hqq = hp.read_map(QQ_file.replace('60.0','60.00'),h=True)
			print('4')
		try:
			uu,huu = hp.read_map(UU_file,h=True)
			print(UU_file)
		except:
			uu,huu = hp.read_map(UU_file.replace('60.0','60.00'),h=True)
			print('5')
		try:
			qu,hqu = hp.read_map(QU_file,h=True)
			print(QU_file)
		except:
			qu,hqu = hp.read_map(QU_file.replace('60.0','60.00'),h=True)
			print('6')
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
				qq_add,hqq_add = hp.read_map(QQ_add.replace('60.0','60.00'),h=True)
			try:
				uu_add,huu_add = hp.read_map(UU_add,h=True)
			except:
				uu_add,huu_add = hp.read_map(UU_add.replace('60.0','60.00'),h=True)
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
	bin_hdu.header['INDXSCHM']='IMPLICIT'
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

# QUIJOTE
mapdir = directory+'quijote_202103_tqu_v1.5/'
noisedir = directory+'quijote_202103_tqu_noise_v1.0/'
noisedir2 = directory+'quijote_202103_sims/'
outdirectory = directory+"quijote_202103_tqu_v1.5_noise_v1.0/"
os.makedirs(outdirectory, exist_ok=True)

for nside in output_nside:
	if nside <= 512:
		namestrings = ['QUIJOTEMFI1_11.0_2021','QUIJOTEMFI1_13.0_2021','QUIJOTEMFI3_11.0_2021','QUIJOTEMFI3_13.0_2021','QUIJOTEMFI2_17.0_2021','QUIJOTEMFI2_19.0_2021','QUIJOTEMFI4_17.0_2021','QUIJOTEMFI4_19.0_2021']
		for i in range(0,len(namestrings)):
			try:
				if i == 0 or i == 1:
					docombine(outdirectory+str(nside)+'_60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits',\
					mapdir+str(nside)+'_60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits',\
					noisedir+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_'+str(nside)+'.fits',comment=comment)
				else:
					docombine(outdirectory+str(nside)+'_60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits',\
					mapdir+str(nside)+'_60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits',\
					noisedir+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_Q_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_U_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_QU_'+str(nside)+'.fits',comment=comment)
					docombine(outdirectory+str(nside)+'_60.0smoothed_'+namestrings[i]+'simnoise_mKCMBunits.fits',\
					mapdir+str(nside)+'_60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits',\
					noisedir2+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_'+str(nside)+'.fits',\
					noisedir2+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_Q_'+str(nside)+'.fits',\
					noisedir2+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_U_'+str(nside)+'.fits',\
					noisedir2+'60.0smoothed_'+namestrings[i]+'_mKCMBunits.fits_variance_QU_'+str(nside)+'.fits',comment=comment)
			except:
				continue

# WMAP9
mapdir = directory+'wmap9_tqu_v1.5/'
noisedir = directory+'wmap9_tqu_noise_v1.0/'
outdirectory = directory+"wmap9_tqu_v1.5_noise_v1.0/"
os.makedirs(outdirectory, exist_ok=True)

for nside in output_nside:
	if nside <= 512:
		namestrings = ['22.8_512_2013','33.0_512_2013','40.7_512_2013','60.7_512_2013','93.5_512_2013']
		for i in range(0,len(namestrings)):
			try:
				docombine(outdirectory+str(nside)+'_60.0smoothed_wmap9beamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
					mapdir+str(nside)+'_60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment)
			except:
				continue

# WMAP9 10k
mapdir = directory+'wmap9_tqu_v1.5/'
noisedir = directory+'wmap9_tqu_noise_v1.0_10k/'
outdirectory = directory+"wmap9_tqu_v1.5_noise_v1.0_10k/"
os.makedirs(outdirectory, exist_ok=True)

for nside in output_nside:
	if nside <= 512:
		namestrings = ['22.8_512_2013','33.0_512_2013','40.7_512_2013','60.7_512_2013','93.5_512_2013']
		for i in range(0,len(namestrings)):
			try:
				docombine(outdirectory+str(nside)+'_60.0smoothed_wmap9beamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
					mapdir+str(nside)+'_60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
					noisedir+'60.0smoothed_wmap9beam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment)
			except:
				continue

# Planck 2015
mapdir = directory+'planck2015_tqu_v1.5/'
noisedir = directory+'planck2015_tqu_noise_v1.0/'
outdirectory = directory+"planck2015_tqu_v1.5_noise_v1.0/"
os.makedirs(outdirectory, exist_ok=True)
rescale=1000.0**2
for nside in output_nside:
	# Standard maps
	namestrings = ['28.4_1024_2015','44.1_1024_2015','70.4_2048_2015','100_2048_2015','143_2048_2015','217_2048_2015','353_2048_2015']
	namestrings2 = ['28.4_1024_2015','44.1_1024_2015','70.4_1024_2015','100_1024_2015','143_1024_2015','217_1024_2015','353_1024_2015']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR2fullbeamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR2fullbeam_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Bandpass subtracted maps
	namestrings = ['28.4_256_2015','44.1_256_2015','70.4_256_2015']
	namestrings2 = ['28.4_1024_2015','44.1_1024_2015','70.4_1024_2015']
	namestrings3 = ['30', '44', '70']
	if nside <= 256:
		for i in range(0,len(namestrings)):
			# try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR2fullbeambpcorrNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR2fullbeambpcorr_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale,	QQ_add = noisedir+'60.0smoothed_PlanckR2bpasscorrection_'+namestrings3[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',	UU_add = noisedir+'60.0smoothed_PlanckR2bpasscorrection_'+namestrings3[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits')
			# except:
			# 	continue

	# I only maps
	namestrings = ['545_2048_2015','857_2048_2015']
	namestrings2 = ['545_1024_2015','857_1024_2015']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR2fullbeamNoise_'+namestrings[i]+'_MJySrunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR2fullbeam_'+namestrings[i]+'_MJySrunits.fits',\
				noisedir+'60.0smoothed_PlanckR2fullbeam_'+namestrings2[i]+'_MJySrunits_variance_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue

# Planck 2018
mapdir = directory+'planck2018_tqu_v1.5/'
noisedir = directory+'planck2018_tqu_noise_v1.0/'
outdirectory = directory+"planck2018_tqu_v1.5_noise_v1.0/"
os.makedirs(outdirectory, exist_ok=True)
for nside in output_nside:
	# Standard maps
	namestrings = ['28.4_1024_2018','44.1_1024_2018','70.4_1024_2018','100_2048_2018','143_2048_2018','217_2048_2018','353_2048_2018']
	namestrings2 = ['28.4_1024_2018','44.1_1024_2018','70.4_1024_2018','100_1024_2018','143_1024_2018','217_1024_2018','353_1024_2018']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR3fullbeamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR3fullbeam_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Intensity only maps
	namestrings = ['545_2048_2018','857_2048_2018']
	namestrings2 = ['545_1024_2018','857_1024_2018']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR3fullbeamNoise_'+namestrings[i]+'_MJySrunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR3fullbeam_'+namestrings[i]+'_MJySrunits.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_MJySrunits_variance_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Bandpass unsubtracted maps
	namestrings = ['28.4_1024_2018','44.1_1024_2018','70.4_1024_2018']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR3fullbeamnobpNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR3fullbeamnobp_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeamnobp_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeamnobp_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeamnobp_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeamnobp_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# deconvolved maps
	namestrings = ['dec40_28.4_1024_2018','dec30_44.1_1024_2018','dec20_70.4_1024_2018']
	namestrings2 = ['28.4_1024_2018','44.1_1024_2018','70.4_1024_2018']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR3decNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR3'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue

mapdir = directory+'planck2018_tqu_v1.5/'
noisedir = directory+'planck2018_tqu_noise_v1.0_10k/'
outdirectory = directory+"planck2018_tqu_v1.5_noise_v1.0_10k/"
os.makedirs(outdirectory, exist_ok=True)
for nside in output_nside:
	# Planck 2018 10k
	# Standard maps
	namestrings = ['28.4_1024_2018','44.1_1024_2018']
	namestrings2 = ['28.4_1024_2018','44.1_1024_2018']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR3fullbeamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR3fullbeam_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings2[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue



# Planck 2020
mapdir = directory+'planck2020_tqu_v1.5/'
noisedir = directory+'planck2020_tqu_noise_v1.0/'
outdirectory = directory+"planck2020_tqu_v1.5_noise_v1.0/"
os.makedirs(outdirectory, exist_ok=True)
for nside in output_nside:
	# Standard maps
	namestrings = ['28.4_1024_2020','44.1_1024_2020','70.4_1024_2020','100_2048_2020','143_2048_2020','217_2048_2020','353_2048_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR4fullbeamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Dipole subtracted maps - using the same noise realisations
	namestrings = ['28.4_1024_2020','44.1_1024_2020','70.4_1024_2020','100_2048_2020','143_2048_2020','217_2048_2020','353_2048_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR4fullbeamnodpNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR4fullbeamnodp_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Intensity only maps
	namestrings = ['545_2048_2020','857_2048_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR4fullbeamNoise_'+namestrings[i]+'_MJySrunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_MJySrunits.fits',\
				noisedir+'60.0smoothed_PlanckR3fullbeam_'+namestrings[i]+'_MJySrunits_variance_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue


# Planck 2020 10k
mapdir = directory+'planck2020_tqu_v1.5/'
noisedir = directory+'planck2020_tqu_noise_v1.0_10k/'
outdirectory = directory+"planck2020_tqu_v1.5_noise_v1.0_10k/"
os.makedirs(outdirectory, exist_ok=True)
for nside in output_nside:
	# Standard maps
	namestrings = ['28.4_1024_2020','44.1_1024_2020','70.4_1024_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR4fullbeamNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
	# Dipole subtracted maps - using the same noise realisations
	namestrings = ['28.4_1024_2020','44.1_1024_2020','70.4_1024_2020']
	for i in range(0,len(namestrings)):
		try:
			docombine(outdirectory+str(nside)+'_60.0smoothed_PlanckR4fullbeamnodpNoise_'+namestrings[i]+'_mKCMBunits.fits',\
				mapdir+str(nside)+'_60.0smoothed_PlanckR4fullbeamnodp_'+namestrings[i]+'_mKCMBunits.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_Q_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_U_'+str(nside)+'.fits',\
				noisedir+'60.0smoothed_PlanckR4fullbeam_'+namestrings[i]+'_mKCMBunits_variance_QU_'+str(nside)+'.fits',comment=comment,rescale=rescale)
		except:
			continue
			continue

