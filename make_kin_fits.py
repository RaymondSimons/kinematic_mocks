'''
Generate IFU maps from SUNRISE output, need mcrx.fits
created: 08/18/2016, Raymond Simons
'''
import astropy
from astropy.io import fits
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from joblib import Parallel, delayed
import glob
from numpy import ma
import sys
from astropy import constants
import astropy
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel, convolve_fft, convolve
import scipy
import os, sys, argparse
from scipy.optimize import curve_fit
import pickle
import pyfits
import os
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot  import *
import numpy as np
from numpy import *


plt.ioff()
def gauss(x, *p):
    A, mu, sigma, A0= p
    if A < 0:
        return nan*zeros(len(x))
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))+A0


def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''\
                                Generate the cameras to use in Sunrise and make projection plots
                                of the data for some of these cameras. Then export the data within
                                the fov to a FITS file in a format that Sunrise understands.
                                ''')
     
 
    parser.add_argument('-a', '--scale', default=400, type=int,
                        help='Scale factor')
    
    args = vars(parser.parse_args())
    return args


class kin_map():
    def __init__(self, data, header, vel_kms, lam,  camera, scale):
        self.orig_cube     = data
        self.orig_cube_hdr = header.copy()
        self.cube          = data.copy()
        self.cube_hdr      = header.copy()
        self.zsize         = data.shape[0]
        self.xsize         = data.shape[1]
        self.ysize         = data.shape[2]
        self.vscale        = vel_kms
        self.lam           = lam
        self.camera        = camera
        self.ascale        = scale

    def rebin(self, new_shape):
        #The original cube has 400 pixels, which is not necessarily evenly divisible by our requested size
        #Necessary size of cube
        temp_orig_shape = new_shape[0]*ceil(self.orig_cube.shape[1]/float(new_shape[0]))
        print temp_orig_shape, new_shape[0], self.orig_cube.shape[1]

        self.cube = zeros((self.zsize, new_shape[0],new_shape[1]))
        self.cube_hdr['CD1_1'] *=  (temp_orig_shape)/new_shape[0]
        self.cube_hdr['CD2_2'] *=  (temp_orig_shape)/new_shape[1]
        self.cube_hdr['CRPIX1'] =  new_shape[0]/2.
        self.cube_hdr['CRPIX2'] =  new_shape[1]/2.
        self.cube_hdr['NAXIS1'] =  new_shape[0]
        self.cube_hdr['NAXIS2'] =  new_shape[1]
        self.cube_hdr['EXTNAME'] =  self.cube_hdr['EXTNAME'].replace('_ORIG', '_OBS')

        pixel_extension = int((temp_orig_shape - self.orig_cube.shape[1])/2.)

        print pixel_extension

        for i in arange(self.zsize):
            M, N = temp_orig_shape, temp_orig_shape
            m, n = new_shape
            temp_orig_cube_slice = zeros((temp_orig_shape, temp_orig_shape))
            x0, y0 = pixel_extension, pixel_extension
            x1, y1 = temp_orig_shape - pixel_extension, temp_orig_shape - pixel_extension
            temp_orig_cube_slice[x0:x1, y0:y1] = self.orig_cube[i]
            if m<M:
                #old way: self.cube[i] = self.orig_cube[i].reshape((m,M/m,n,N/n)).mean(3).mean(1)
                #new way:
                self.cube[i] = temp_orig_cube_slice.reshape((m,M/m,n,N/n)).mean(3).mean(1)
            else:
                #old way: self.cube[i] = np.repeat(np.repeat(self.orig_cube[i], m/M, axis=0), n/N, axis=1)
                self.cube[i] = np.repeat(np.repeat(temp_orig_cube_slice, m/M, axis=0), n/N, axis=1)



        self.zsize      = self.cube.shape[0]
        self.xsize      = self.cube.shape[1]
        self.ysize      = self.cube.shape[2]

    def generate_blurred_map(self, kernel_size_arc = 0.6, band = 'H'):
        self.blrcube    = self.cube.copy()*nan

        #The pixel values in our cube are W/m/m^2/Sr (surface brightness). To get to units of 
        #flux density per pixel
        #check to make sure this pixel scale is in proper coordinates
        pix_scale_kpc = self.cube_hdr['CD1_1']*u.kpc/u.pixel
        print pix_scale_kpc #kpc per pixel
        pix_scale_arc = pix_scale_kpc * cosmo.arcsec_per_kpc_proper(1./self.ascale-1)
        print pix_scale_arc #arc per pixel
        pix_scale_str = (pix_scale_arc**2.).to(u.steradian/u.pixel**2.)
        print pix_scale_str #steradian per square pixel


        #set kernel size in pixels
        self.kernel_size_arc = kernel_size_arc*u.arcsec
        self.kernel_size_pix = self.kernel_size_arc/pix_scale_arc


        print 'Seeing:'
        print '\t sigma = ', self.kernel_size_pix, self.kernel_size_arc 
        print '\t fwhm = ', 2.35*self.kernel_size_pix, 2.35 * kernel_size_arc

        #Generate the kernel from the seeing size in pixels
        self.kernel = Gaussian2DKernel(self.kernel_size_pix.value)

        print 'Convolving spatially...'
        for i in arange(self.zsize):
            self.blrcube[i] = convolve_fft(self.cube[i], self.kernel)
        #self.blrcube += np.random.normal(0, max(self.blrcube.ravel())/100., shape(self.blrcube))

        #KMOS reaches a point source 5-sigma sensitvity in 8 hr of
        #of (J, H, K) = (22, 21.0, 20.5) AB magnitudes 
        #for R ~ (3380, 3800, 3750)
        #baseline sensitivity measurements from: http://www2011.mpe.mpg.de/Highlights/FB2004/exp13_bender.pdf
        if   band == 'H': sens, R = 21.0, 3800
        elif band == 'J': sens, R = 22.0, 3380
        elif band == 'K': sens, R = 20.5, 3750
        else: 
            sens, R = 21.0, 3800.
            print 'Bad input band, setting sensitivty to H-band value, %s AB mag'%(sens)


        print 'Convolving spectrally...'

        print 'Adding noise...'

        sens_si_fd = (sens*u.ABmag).to(u.Watt/(u.meter*u.meter)/(u.Hz))*astropy.constants.c/(self.lam[0]*u.meter)**2.
        print sens_si_fd


        #Currently in m, want to get in terms of hz^-1. F_v = (F_lam)*lam^2/c
        #Ths units of this factor are 1/(str*m). Let's now consider an 'aperture' equal to 
        #1 spatial fwhm * 1 spectral fwhm. We want to add noise equal to 1/5th the sensitivity (i.e., 1 sigma sensitivty) over this aperture.
        #In steradians, the PSF is:

        psf_str = pi*(((2.35/2.)*self.kernel_size_arc)**2.).to(u.steradian/u.pixel**2.)
        print '\t', psf_str, 'steradian seeing FWHM area'


        #The spectral lsf fwhm (in pixels) is:
        self.kms_per_pix = self.vscale[1]-self.vscale[0]
        lsf_pix = (3.e5/R)/self.kms_per_pix
        print lsf_pix
        
        sens_noise = sens_si_fd/psf_str/lsf_pix
        print sens_noise
        self.blrcube += np.random.normal(0, sens_noise.value, self.cube.shape)

    def generate_intrinsic_kin_map(self):
        self.int_vel_hdr = self.orig_cube_hdr.copy()
        self.int_vel_hdr['IMUNIT'] = 'km/s'
        self.int_vel_hdr.remove('UNITCONV')


        self.vel_int    = np.zeros((self.xsize, self.ysize))*np.nan
        self.evel_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.disp_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.edisp_int  = np.zeros((self.xsize, self.ysize))*np.nan
        self.ha_int     = np.zeros((self.xsize, self.ysize))*np.nan

        cb_std = self.orig_cube[0:5,:,:].std()
        
        #shutting off for testing
        '''
        for i in arange(self.cube.shape[1]):
            for j in arange(self.cube.shape[2]):  
                if self.cube[:,i,j].max() > 10*cb_std: #likely have an Ha line here
                    spec = self.cube[:,i,j]
                    try:
                        c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), self.vscale[nanargmax(spec)], 10, spec[0]])  
                        #if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                        #(c_a[0] > 0) & (c_a[2] > 5) & (sqrt(v_a[2,2]) < 30) & (sqrt(v_a[1,1]) < 30):
                        self.vel_int[i,j]   = c_a[1]
                        self.evel_int[i,j]  = sqrt(v_a[1,1])
                        self.disp_int[i,j]  = c_a[2]
                        self.edisp_int[i,j] = sqrt(v_a[2,2])
                        self.ha_int[i,j] = c_a[0]*c_a[2]*sqrt(2*pi)                            
                    except:
                        pass
        '''
    def generate_observed_kin_map(self):
        self.obs_vel_hdr = self.cube_hdr.copy()
        self.obs_vel_hdr['IMUNIT'] = 'km/s'
        self.obs_vel_hdr.remove('UNITCONV')

        self.vel_obs    = np.zeros((self.xsize, self.ysize))*np.nan
        self.evel_obs   = np.zeros((self.xsize, self.ysize))*np.nan
        self.disp_obs   = np.zeros((self.xsize, self.ysize))*np.nan
        self.edisp_obs  = np.zeros((self.xsize, self.ysize))*np.nan
        self.ha_obs     = np.zeros((self.xsize, self.ysize))*np.nan
        for i in arange(self.blrcube.shape[1]):
            for j in arange(self.blrcube.shape[2]):                
                krnl = Gaussian1DKernel(1)
                spec = convolve_fft(self.blrcube[:,i,j], krnl)
                try:
                    c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), self.vscale[nanargmax(spec)], 30, median(spec[0:15])])           
                    #if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                    #(c_a[0] > 0) & (c_a[2] > 10) & (sqrt(v_a[2,2]) < 30) & (sqrt(v_a[1,1]) < 30):
                    self.vel_obs[i,j]   = c_a[1]
                    self.evel_obs[i,j]  = sqrt(v_a[1,1])
                    self.disp_obs[i,j]  = c_a[2]
                    self.edisp_obs[i,j] = sqrt(v_a[2,2])
                    self.ha_obs[i,j] = c_a[0]*c_a[2]*sqrt(2*pi)                            
                except:
                    pass
        #self.vel_obs[-isnan(self.vel_obs)] += np.random.normal(0,5,shape(self.vel_obs[-isnan(self.vel_obs)]))
        #self.disp_obs[-isnan(self.disp_obs)] += np.random.normal(0,5,shape(self.vel_obs[-isnan(self.vel_obs)]))

    def get_hdulist(self, master_hdulist):
        colhdr = fits.Header()
        #master_hdulist.append(fits.ImageHDU(data = self.orig_cube, header = self.orig_cube_hdr, name = 'cam%i_orig_cube'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.cube, header = self.cube_hdr, name = 'cam%i_cub_int'%self.camera))        
        master_hdulist.append(fits.ImageHDU(data = self.blrcube, header = self.cube_hdr, name = 'cam%i_cub_obs'%self.camera))

        master_hdulist.append(fits.ImageHDU(data = array([self.disp_int, self.edisp_int]),header = self.int_vel_hdr,  name = 'cam%i_dis_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_obs, self.edisp_obs]), header = self.obs_vel_hdr,name = 'cam%i_dis_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_int,self.evel_int]),header = self.int_vel_hdr,  name = 'cam%i_vel_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_obs,self.evel_obs]), header = self.obs_vel_hdr, name = 'cam%i_vel_obs'%self.camera))

        self.ha_int_hdr = self.orig_cube_hdr
        self.ha_int_hdr.remove('IMUNIT')
        self.ha_int_hdr.remove('UNITCONV')


        self.ha_obs_hdr = self.cube_hdr
        self.ha_obs_hdr.remove('IMUNIT')
        self.ha_obs_hdr.remove('UNITCONV')


        master_hdulist.append(fits.ImageHDU(data = self.ha_int, header = self.orig_cube_hdr, name = 'cam%i_hal_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_obs, header = self.cube_hdr, name = 'cam%i_hal_obs'%self.camera))

        return master_hdulist


def run_kin_fits(abspath, scale, kmap_name, gal, outdir, mcrx_data):
    print '\tReading in mcrx file for (%s, %.3f)'%(gal, scale)

    #setting constants
    Ha_m = 6.563e-7 #Halpha in meters
    c_kms = constants.c.value*1.e-3    #speed of light in km/s

    #reading in data
    #mcrx_data = fits.open(abspath)

    #reading number of cameras
    ncams = mcrx_data['MCRX'].header['N_CAMERA']


    #generate empty fits file
    master_hdulist = []

    #generate primary header for fits file
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Storing the kinematic maps in this FITS file."
    prihdr['ncams'] = str(ncams)
    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)


    #reading in lambda and velocity arrays
    lam = mcrx_data['LAMBDA'].data['lambda']
    vel_pixel = c_kms*(lam[1]-lam[0])/lam[0]
    vel_arr = vel_pixel * arange(len(lam))
    loc_Ha = argmin(abs(lam-Ha_m))
    vel_arr-=vel_arr[loc_Ha]

    #add the velocity and lambda arrays in the first header
    cols = fits.ColDefs([fits.Column(name='velocity', format = 'D',array=vel_arr, unit = 'km/s'), 
                         fits.Column(name='lambda', format = 'D', array=lam, unit = 'm')])
    colhdr = fits.Header()
    master_hdulist.append(fits.BinTableHDU.from_columns(cols, name = 'props', header = colhdr))


    #run kinematic fitting routine for all cameras
    #for cam_n in arange(ncams):
    for cam_n in arange(10,11):
        print '\t\t Running on (%s, %.3f, %i)'%(gal, scale, cam_n)
        camera = mcrx_data['CAMERA%i'%(cam_n)]   
        kmap = kin_map(camera.data, camera.header, vel_arr, lam,  cam_n, scale)
        print 'Generating intrinsic kinematic map'
        kmap.generate_intrinsic_kin_map()

        kmap.rebin([60,60])
        kmap.generate_blurred_map(kernel_size_arc = 0.6)
        kmap.generate_observed_kin_map()
        master_hdulist = kmap.get_hdulist(master_hdulist)

    #Save the fits file
    thdulist = fits.HDUList(master_hdulist)
    print '\tSaving to ' + outdir+'/'+kmap_name
    thdulist.writeto(outdir+'/'+kmap_name, clobber = True)


if __name__ == '__main__':
    path_to_mcrx = './*/ifu/'
    mcrx_files = glob.glob(path_to_mcrx+'/mcrx.fits.gz')
    scales = []
    abspaths = []
    kmap_names = []
    gal = os.path.abspath('.').split('/')[-1]

    print 'Working on %s'%gal


    for n, fl in enumerate(mcrx_files):
        abspaths.append(os.path.abspath(fl))
        sc_loc = abspaths[n].find('_a')
        scales.append(float(abspaths[n][sc_loc+2:sc_loc+7]))
        kmap_names.append('%s_a%.3f_kmap.fits'%(gal, scales[n]))

    #Where to write the kinematic map files
    outdir = '/nobackupp2/rcsimons/data/kin_maps/%s'%gal

    test = True
    if test:
        #want to select individual systems
        scales   = array(scales)
        n_sel = where(scales == 0.38)[0][0]
        #mcrx_data = fits.open(abspaths[n_sel])
        run_kin_fits(abspaths[n_sel], scales[n_sel], kmap_names[n_sel], gal, outdir, mcrx_data)
    else:
        #run on all
        Parallel(n_jobs = -1)(delayed(run_kin_fits)(abspaths[i], scales[i], kmap_names[i], gal, outdir) for i in arange(len(scales)))
    


























