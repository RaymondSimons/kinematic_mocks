import pyfits
import glob
from numpy import *
import os


files = ['VELA28_a0.%3i_sunrise'%(num) for num in arange(300,380,10)]
simdir = '/nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/VELA28'
for fname in files:
	fname_fits = fname+'_cam1.fits'
	print simdir+'/'+fname_fits+'/ifu/mcrx.fits'
	fle = glob.glob(simdir+'/'+fname+'/ifu/mcrx.fits')[0]
	
	data = mcrx['CAMERA1-NONSCATTER'].data
	
	pyfits.writeto('/nobackupp2/rcsimons/sunrise_testing/jwst/'+fname_fits, data, clobber = True)
	#os.system('cp '+fle+' /nobackupp2/rcsimons/sunrise_testing/jwst/'+fname+'_mcrx.fits')

	print '/nobackupp2/rcsimons/sunrise_testing/jwst/'+fname_fits
	'''
	'''