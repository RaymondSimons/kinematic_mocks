import pyfits
import glob
from numpy import *

files = ['VELA28_a0.%3i_sunrise'%(num) for num in arange(300,380,10)]
simdir = '/nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/VELA28'
for fname in files:
	fname_fits = fname+'_cam1.fits'
	print simdir+'/'+fname_fits+'/ifu/mcrx.fits'
	fle = glob.glob(simdir+'/'+fname_fits+'/ifu/mcrx.fits')
	print fle
	mcrx = pyfits.open(fle)
	data = mcrx['CAMERA0-NONSCATTER'].data[1]

	pyfits.writeto('/nobackupp2/rcsimons/sunrise_testing/jwst', data, clobber = True)


	print '/nobackupp2/rcsimons/sunrise_testing/jwst/'+fname_fits
	'''
	'''