import pyfits
import glob
import astropy
from astropy import units
from astropy.cosmology import Planck15 as cosmo
from astropy.convolution import Gaussian2DKernel, convolve_fft, Gaussian1DKernel
plt.close('all')
plt.ioff()

files = glob.glob('../data/jwst/*340*cam1.fits')


a = array([float(file.rsplit('_')[1].replace('a','')) for file in files])
z = 1./a - 1.

arc_per_kpc = cosmo.arcsec_per_kpc_proper(z)


"""
Nirspec design
0.2 x 0.5" shutters
R ~2700

"""
ysize = 1.5
xsize = 0.2

msa_arc_pix = 0.065
mos_arc_pix = 0.18
msa_fwhm_arc = 0.1
mos_fwhm_arc = 0.7

def rebin(a, new_shape):
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)


fig = figure()
fig2 = figure(figsize = (10,5))

ax1 = fig.add_subplot(122)
ax2 = fig.add_subplot(121)

ax3 = fig2.add_subplot(212)
ax4 = fig2.add_subplot(211)

mos_sens = 1.e-17/4./4.1067600250214074e-11
msa_sens = 5.e-19/4./4.1067600250214074e-11



for i, file in enumerate(files):
	data = pyfits.open(file)[0].data[225:325]
	arc_per_pixel = 50./shape(data)[1] * arc_per_kpc[i].value
	data = data * units.W.in_units('erg/s')*1.e-4/units.steradian.in_units('arcsec**2.')

	data_msa = zeros((shape(data)[0], 100,100))
	data_mos = zeros((shape(data)[0], 40,40))

	gauss_kernel_msa = Gaussian2DKernel(msa_fwhm_arc/2.35/msa_arc_pix)
	gauss_kernel_mos = Gaussian2DKernel(mos_fwhm_arc/2.35/mos_arc_pix)

	for j in arange(len(data)):
		data_msa[j] = rebin(data[j], (100,100))
		data_msa[j] = convolve_fft(data_msa[j], gauss_kernel_msa)
		data_mos[j] = rebin(data[j], (40,40))
		data_mos[j] = convolve_fft(data_mos[j], gauss_kernel_mos)

	arc_pix_mos = arc_per_pixel*10
	arc_pix_msa = arc_per_pixel*4

	mos_sens = 1.e-17/4./4.1067600250214074e-11
	msa_sens = 5.4e-19/2./4.1067600250214074e-11*5.



	ymax_msa = shape(data_msa)[1]/2.+xsize/msa_arc_pix/2.
	ymin_msa = shape(data_msa)[1]/2.-xsize/msa_arc_pix/2.
	
	ymax_mos = shape(data_mos)[1]/2.+xsize/mos_arc_pix/2.
	ymin_mos = shape(data_mos)[1]/2.-xsize/mos_arc_pix/2.


	ymax_orig = shape(data)[1]/2.+xsize/arc_per_pixel/2.
	ymin_orig = shape(data)[1]/2.-xsize/arc_per_pixel/2.


	slit_data_msa = sum(data_msa[:,ymin_msa:ymax_msa,:], axis = 1)*arc_pix_msa**2.
	msa_lsf = Gaussian1DKernel(3)

	for q in arange(shape(slit_data_msa)[0]):
		slit_data_msa[q] = convolve_fft(slit_data_msa[q], msa_lsf)


	slit_data_msa+= np.random.normal(0,msa_sens,shape(slit_data_msa))


	ax3.imshow(slit_data_msa.transpose(), cmap = 'Greys_r', interpolation = 'nearest')
	ax3.set_ylim(45-1.25/msa_arc_pix, 45+1.25/msa_arc_pix)
	#ax3.set_xticklabels([''])
	#ax3.set_yticklabels([''])

	slit_data_mos = sum(data_mos[:,ymin_mos:ymax_mos,:], axis = 1)*arc_pix_mos**2.
	mos_lsf = Gaussian1DKernel(3)
	for q in arange(shape(slit_data_msa)[0]):
		slit_data_mos[q] = convolve_fft(slit_data_mos[q], mos_lsf)

	slit_data_mos+= np.random.normal(0,mos_sens,shape(slit_data_mos))
	ax4.imshow(slit_data_mos.transpose(), cmap = 'Greys_r', interpolation = 'nearest')
	ax4.set_ylim(18-1.50/mos_arc_pix, 18+1.50/mos_arc_pix)
	ax4.set_xlim(19,73)
	#ax4.set_xticklabels([''])
	ax3.set_xticks([15,45,75])
	ax3.set_xlim(10,80)
	ax3.set_xticklabels(['1.926', '1.930', '1.934'])
	ax3.set_yticklabels([''])
	ax4.set_yticklabels([''])
	ax4.set_xticklabels([''])
	ax3.set_xlabel(r'$\lambda (\mu m$)', fontsize = 18)
	ax3.set_ylabel('position along slit', fontsize = 14)
	#plt.axis('equal')






















	data = pyfits.open(file)[0].data[225:325]

	ha_map = sum(data, axis = 0)
	ha_map[ha_map == 0] = 0.01
	ax1.imshow(log10(ha_map), cmap = 'Greys', vmin = 0, vmax = 5)
	ax1.axhline(ymax_orig, color = 'red')
	ax1.axhline(ymin_orig, color = 'red')

	for i in arange(0,400,0.5/arc_per_pixel):
		ax1.axvline(i, ymin = ymax_orig/shape(data)[1], ymax = ymin_orig/shape(data)[1], color = 'red')

	ax1.set_xlim(100,280)
	ax1.set_ylim(110,290)
	ax2.imshow(sum(data_new, axis = 0), cmap = 'Greys_r', vmin = 0, vmax = 1000)
	ax2.axhline(ymax)
	ax2.axhline(ymin)


	ax1.set_xticklabels([''])
	ax1.set_yticklabels([''])
	ax3.tick_params(which='major', axis = 'x', width = 2, length=4, color = 'white')
	ax1.tick_params(which='major', axis = 'both', width = 2, length=0, color = 'white')
	ax1.tick_params(which='minor', axis = 'both', width = 2, length=0, color = 'white')



fig2.subplots_adjust(hspace =0.05, bottom = 0.12)
fig2.savefig('../figures/jwst.png', dpi = 200)
fig.savefig('../figures/jwst_ha.png', dpi = 300)

'''
	imshow(sum(data.data, axis = 0), cmap = 'Greys_r', vmin = 0, vmax = 1000)
	slitdata = sum(data.data[200:400, 175:225, 175:225], axis = 1)
	cont = median(slitdata[:,25])
	back = std(slitdata[:,0:10].ravel())

	slit_correct = slitdata-bckrd+np.random.normal(0,400,shape(slitdata))
	#hist((slit_correct).ravel(), bins = linspace(-20,20,100))
	#imshow(slit_correct, cmap='Greys_r', vmin = -3, vmax = 3)
	#imshow(slitdata- median(slitdata), vmin = 0, vmax = 10, cmap = 'Greys_r')
	raw_input('')
	clf()
	plt.close('all')
'''






