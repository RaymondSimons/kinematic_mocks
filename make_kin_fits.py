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
        self.redshift      = 1./scale - 1
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

    def rebin_and_dim(self, new_shape):
        #The original cube has 400 pixels, which is not necessarily evenly divisible by our requested size
        #Necessary size of cube
        print '\t\t Rebinning to ', new_shape 
        temp_orig_shape = new_shape[0]*ceil(self.orig_cube.shape[1]/float(new_shape[0])/2.)*2.
        print temp_orig_shape
        #Get to an even number
        #temp_orig_shape = math.ceil(temp_orig_shape / 2.) * 2
        print temp_orig_shape/float(new_shape[0])


        #print temp_orig_shape, new_shape[0], self.orig_cube.shape[1]

        self.cube = zeros((self.zsize, new_shape[0],new_shape[1]))
        
        self.cube_hdr['CD1_1'] *=  (temp_orig_shape)/new_shape[0]
        self.cube_hdr['CD2_2'] *=  (temp_orig_shape)/new_shape[1]
        self.cube_hdr['CRPIX1'] =  new_shape[0]/2.
        self.cube_hdr['CRPIX2'] =  new_shape[1]/2.
        self.cube_hdr['NAXIS1'] =  new_shape[0]
        self.cube_hdr['NAXIS2'] =  new_shape[1]
        self.cube_hdr['EXTNAME'] =  self.cube_hdr['EXTNAME'].replace('_ORIG', '_OBS')

        pixel_extension = int((temp_orig_shape - self.orig_cube.shape[1])/2.)

        #print pixel_extension

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

        #Surface brightness dimming
        #print 'Applying cosmological surface brightness dimming', already have a term of ^2 from pixel shift
        self.cube = self.cube/(self.redshift + 1.)**5. #apparently it's five for la_nu
        #self.cube = self.cube/(self.redshift + 1.)**4.

        print (self.redshift + 1.)**2.


        self.cube_hdr['orig_linear_fov'] = (self.cube_hdr['linear_fov'], '[kpc] original sunrise linear field of view')
        self.cube_hdr['linear_fov']      = (self.ysize * self.cube_hdr['CD1_1'], '[kpc] Linear field of view in y-dir after rebin')

        self.pix_scale_kpc = self.cube_hdr['CD1_1']*u.kpc/u.pixel
        print self.pix_scale_kpc #kpc per pixel
        self.pix_scale_arc = self.pix_scale_kpc * cosmo.arcsec_per_kpc_proper(1./self.ascale-1)
        self.cube_hdr['pix_size']=(self.pix_scale_arc.value, '[arcsec] per side')

        print self.pix_scale_arc #arc per pixel
        self.pix_scale_str = (self.pix_scale_arc**2.).to(u.steradian/u.pixel**2.)
        print self.pix_scale_str #steradian per square pixel


    def generate_blurred_map(self, kernel_size_arc = 0.6, band = 'H'):
        #KMOS reaches a point source 5-sigma sensitvity in 8 hr of
        #of (J, H, K) = (22, 21.0, 20.5) AB magnitudes 
        #for R ~ (3380, 3800, 3750)
        #baseline sensitivity measurements from: http://www2011.mpe.mpg.de/Highlights/FB2004/exp13_bender.pdf
        if   band == 'H': sens, R = 24.0, 3800 #sens, R = 21.0, 3800
        if   band == 'H': sens, R = 28.0, 2700 #sens, R = 21.0, 3800
        #if   band == 'H': sens, R = 26.0, 3800 #sens, R = 21.0, 3800
        elif band == 'J': sens, R = 22.0, 3380
        elif band == 'K': sens, R = 20.5, 3750
        else: 
            sens, R = 21.0, 3800.
            print 'Bad input band, setting sensitivty to H-band value, %s AB mag'%(sens)

        self.blrcube    = self.cube.copy()*nan

        #The pixel values in our cube are W/m/m^2/Sr (surface brightness). 
        


        #Generate the kernel from the seeing size in pixels

        print 'Convolving spatially...'
        if True:
            #set kernel size in pixels
            self.kernel_size_arc = kernel_size_arc*u.arcsec
            self.kernel_size_pix = self.kernel_size_arc/self.pix_scale_arc
            self.psf_str = pi*(((2.35/2.)*self.kernel_size_arc)**2.).to(u.steradian)
            print '\t\tSeeing:'
            print '\t\t\t sigma = ', self.kernel_size_pix, ',', self.kernel_size_arc 
            print '\t\t\t fwhm = ', 2.35*self.kernel_size_pix, ',', 2.35 * self.kernel_size_arc
            print '\t\t\t fwhm area = ', self.psf_str #in steradians

            self.cube_hdr['seeing']=(self.kernel_size_arc.value, '[%s] fwhm'%str(self.kernel_size_arc.unit))

            self.kernel = Gaussian2DKernel(self.kernel_size_pix.value)
            for i in arange(self.zsize):
                self.blrcube[i] = convolve_fft(self.cube[i], self.kernel)


        print 'Convolving spectrally...'
        if True:
            #The spectral lsf fwhm (in pixels) is:
            self.kms_per_pix = (self.vscale[1]-self.vscale[0])*u.km/u.s/u.pixel
            self.lsf_kms = astropy.constants.c.to('km/s')/R/2.35
            self.lsf_pix = self.lsf_kms/self.kms_per_pix

            print '\t\t Line spread function:'
            print '\t\t\t sigma = ', self.lsf_pix, ',', self.lsf_kms
            print '\t\t\t fwhm = ', 2.35*self.lsf_pix, ',', 2.35 * self.lsf_kms
            self.cube_hdr['LSF']=(2.35*self.lsf_kms.value, '[%s] fwhm'%str(self.lsf_kms.unit))
            self.cube_hdr['R']=(R, 'spectral resolution')

            self.spec_kernel = Gaussian1DKernel(self.lsf_pix.value)

            for xx in arange(self.xsize):
                for yy in arange(self.ysize):
                    self.blrcube[:, xx, yy] = convolve_fft(self.blrcube[:,xx, yy], self.spec_kernel)

        print 'Adding noise...'
        self.cube_hdr['sens']=(sens, '[AB mag] 5 sigma 8 hr sensitivty')

        sens_si_fd = (sens*u.ABmag).to(u.Watt/(u.meter*u.meter)/(u.Hz))*astropy.constants.c/(self.lam[0]*u.meter)**2.
        print '\t\t', sens_si_fd


        #Currently in m, want to get in terms of hz^-1. F_v = (F_lam)*lam^2/c
        #Ths units of this factor are 1/(str*m). Let's now consider an 'aperture' equal to 
        #1 spatial fwhm * 1 spectral fwhm. We want to add noise equal to 1/5th the sensitivity (i.e., 1 sigma sensitivty) over this aperture.
        #In steradians, the PSF is:

        sens_noise = sens_si_fd/self.psf_str/self.lsf_pix/5.
        print '\t\t', sens_noise
        self.cube_hdr['noise'] = (sens_noise.value, '[%s]'%sens_noise.unit)
        self.blrcube += np.random.normal(0, sens_noise.value, self.cube.shape)

    def generate_intrinsic_kin_map(self):
        self.int_vel_hdr = self.orig_cube_hdr.copy()
        self.int_vel_hdr['IMUNIT'] = ('km/s','unit')
        self.int_vel_hdr.remove('UNITCONV')


        self.vel_int    = np.zeros((self.xsize, self.ysize))*np.nan
        self.evel_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.disp_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.edisp_int  = np.zeros((self.xsize, self.ysize))*np.nan
        self.ha_int     = np.zeros((self.xsize, self.ysize))*np.nan

        cb_std = self.orig_cube[0:5,:,:].std()
        
        for i in arange(self.cube.shape[1]):
            for j in arange(self.cube.shape[2]):  
                if self.cube[:,i,j].max() > 5*cb_std: #likely have an Ha line here
                    spec = self.cube[:,i,j]
                    try:
                        c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), self.vscale[nanargmax(spec)], 10, spec[0]])  
                        if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                        (c_a[0] > 0) & (c_a[2] > 0) & (sqrt(v_a[2,2]) < 30) & (sqrt(v_a[1,1]) < 30):
                            self.vel_int[i,j]   = c_a[1]
                            self.evel_int[i,j]  = sqrt(v_a[1,1])
                            self.disp_int[i,j]  = c_a[2]
                            self.edisp_int[i,j] = sqrt(v_a[2,2])
                            self.ha_int[i,j] = c_a[0]*c_a[2]*sqrt(2*pi)                            
                    except:
                        print 'Intrinsic kinematic fit broke at pixel %i %i'%(i,j)
                        pass

    def generate_observed_kin_map(self):
        self.obs_vel_hdr = self.cube_hdr.copy()
        self.obs_vel_hdr['IMUNIT'] = ('km/s', 'unit')
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
                    if (isfinite(sqrt(v_a[1,1]))) & (c_a[2]**2. > self.lsf_kms.value**2.) & (isfinite(sqrt(v_a[2,2]))) & (c_a[0] > 0):
                        self.vel_obs[i,j]   = c_a[1]
                        self.evel_obs[i,j]  = sqrt(v_a[1,1])
                        self.disp_obs[i,j]  = sqrt(c_a[2]**2. - self.lsf_kms.value**2.)
                        self.edisp_obs[i,j] = sqrt(v_a[2,2])
                        self.ha_obs[i,j] = c_a[0]*c_a[2]*sqrt(2*pi)                            
                except:
                    print 'Observed kinematic fit broke at pixel %i %i'%(i,j)
                    pass
        #self.vel_obs[-isnan(self.vel_obs)] += np.random.normal(0,5,shape(self.vel_obs[-isnan(self.vel_obs)]))
        #self.disp_obs[-isnan(self.disp_obs)] += np.random.normal(0,5,shape(self.vel_obs[-isnan(self.vel_obs)]))

    def get_hdulist(self, master_hdulist):
        colhdr = fits.Header()
        #master_hdulist.append(fits.ImageHDU(data = self.orig_cube, header = self.orig_cube_hdr, name = 'cam%i_orig_cube'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.cube, header = self.cube_hdr, name = 'cam%.2i_cub_int'%self.camera))        
        master_hdulist.append(fits.ImageHDU(data = self.blrcube, header = self.cube_hdr, name = 'cam%.2i_cub_obs'%self.camera))

        master_hdulist.append(fits.ImageHDU(data = array([self.disp_int, self.edisp_int]),header = self.int_vel_hdr,  name = 'cam%.2i_dis_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_obs, self.edisp_obs]), header = self.obs_vel_hdr,name = 'cam%.2i_dis_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_int,self.evel_int]),header = self.int_vel_hdr,  name = 'cam%.2i_vel_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_obs,self.evel_obs]), header = self.obs_vel_hdr, name = 'cam%.2i_vel_obs'%self.camera))

        self.ha_int_hdr = self.orig_cube_hdr.copy()
        self.ha_int_hdr['IMUNIT'] = (self.orig_cube_hdr['IMUNIT']+ ' km/s', 'A*sigma*sqrt(2*pi) of gauss fit')
        self.ha_int_hdr.remove('UNITCONV')


        self.ha_obs_hdr = self.cube_hdr.copy()
        self.ha_obs_hdr['IMUNIT'] = (self.orig_cube_hdr['IMUNIT']+ ' km/s', 'A*sigma*sqrt(2*pi) of gauss fit')
        self.ha_obs_hdr.remove('UNITCONV')
        self.ha_obs_hdr.remove('LSF')
        self.ha_obs_hdr.remove('R')


        master_hdulist.append(fits.ImageHDU(data = self.ha_int, header = self.ha_int_hdr, name = 'cam%i_hal_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_obs, header = self.ha_obs_hdr, name = 'cam%i_hal_obs'%self.camera))

        return master_hdulist


def run_kin_fits(abspath, scale, kmap_name, gal, outdir, mcrx_data, arc_per_pixel = 0.1):# testing
#def run_kin_fits(abspath, scale, kmap_name, gal, outdir, arc_per_pixel = 0.2):
    print '\tReading in mcrx file for (%s, %.3f)'%(gal, scale)

    #setting constants
    Ha_m = 6.563e-7 #Halpha in meters
    c_kms = constants.c.value*1.e-3    #speed of light in km/s

    #reading in data
    #mcrx_data = fits.open(abspath) #testing

    #reading number of cameras
    ncams = mcrx_data['MCRX'].header['N_CAMERA']


    #generate empty fits file
    master_hdulist = []

    #generate primary header for fits file
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Storing the kinematic maps in this FITS file."
    prihdr['name'] = gal
    prihdr['filename'] = kmap_name
    prihdr['ncams'] = ncams

    prihdr['z'] = (round(1./scale - 1,4), 'redshift (w/ %s)'%cosmo.name)
    prihdr['ascale'] = (scale, 'scale factor')

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
    for cam_n in arange(10,11): #testing
    #for cam_n in arange(ncams):
        np.random.seed()
        print '\t\t Running on (%s, %.3f, %i)'%(gal, scale, cam_n)
        camera = mcrx_data['CAMERA%i'%(cam_n)]   
        camera_params =  mcrx_data['CAMERA%i-PARAMETERS'%(cam_n)].header
        if True:
            camera.header['z'] = (round(1./scale - 1,4), 'redshift (w/ %s)'%cosmo.name)
            camera.header['ascale'] = (scale, 'scale factor')
            camera.header['camera'] = (cam_n, 'camera')


            camera.header['cameradist'] = (camera_params['cameradist'], '[kpc] Distance from origin to camera')
            camera.header['theta'] = (camera_params['theta'], '[rad] Angular coordinate theta of camera positi')
            camera.header['phi']   = (camera_params['phi'], '[rad] Angular coordinate phi of camera position')

            camera.header['CAMPOSX'] = (camera_params['CAMPOSX'], '[kpc] X position of camera        ')
            camera.header['CAMPOSY'] = (camera_params['CAMPOSY'], '[kpc] Y position of camera        ')
            camera.header['CAMPOSZ'] = (camera_params['CAMPOSZ'], '[kpc] Z position of camera        ')
            camera.header['CAMDIRX'] = (camera_params['CAMDIRX'], '[kpc] X comp of camera viewing dir')
            camera.header['CAMDIRY'] = (camera_params['CAMDIRY'], '[kpc] Y comp of camera viewing dir')
            camera.header['CAMDIRZ'] = (camera_params['CAMDIRZ'], '[kpc] Z comp of camera viewing dir')
            camera.header['CAMUPX']  = (camera_params['CAMUPX'] , '[kpc] X comp of camera Y axis dir ')
            camera.header['CAMUPY']  = (camera_params['CAMUPY'] , '[kpc] Y comp of camera Y axis dir ')
            camera.header['CAMUPZ']  = (camera_params['CAMUPZ'] , '[kpc] Z comp of camera Y axis dir ')
            camera.header['linear_fov'] = (camera_params['linear_fov'], '[kpc] Linear field of view in y-dir at origin') 


        kmap = kin_map(camera.data, camera.header, vel_arr, lam,  cam_n, scale)
        print 'Generating intrinsic kinematic map'
        kmap.generate_intrinsic_kin_map()

        print 'linear fov', kmap.cube_hdr['linear_fov']
        npix_new = ceil((kmap.cube_hdr['linear_fov']*cosmo.arcsec_per_kpc_proper(2).value)/arc_per_pixel/2.)*2.
        kmap.rebin_and_dim([npix_new, npix_new])
        kmap.generate_blurred_map(kernel_size_arc = 0.05)
        kmap.generate_observed_kin_map()
        master_hdulist = kmap.get_hdulist(master_hdulist)

    #Save the fits file
    thdulist = fits.HDUList(master_hdulist)
    print '\tSaving to ' + outdir+'/'+kmap_name
    thdulist.writeto(outdir+'/jwst_'+kmap_name, clobber = True)


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
    if test: #testing
        scales   = array(scales)
        n_sel = where(scales == 0.390)[0][0] #want to select individual systems
        print 'Reading in mcrx file...'
        #mcrx_data = fits.open(abspaths[n_sel])
        run_kin_fits(abspaths[n_sel], scales[n_sel], kmap_names[n_sel], gal, outdir, mcrx_data)
    else:
        #run on all
        Parallel(n_jobs = -1)(delayed(run_kin_fits)(abspaths[i], scales[i], kmap_names[i], gal, outdir) for i in arange(len(scales)))
        #Parallel(n_jobs = -1)(delayed(run_kin_fits)(abspaths[i], scales[i], kmap_names[i], gal, outdir) for i in arange(10, 15,1)) #testing
    


























