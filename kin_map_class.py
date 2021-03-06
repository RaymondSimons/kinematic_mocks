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


class kin_map():
    def __init__(self, data, header, vel_kms, lam,  camera, scale):
        self.orig_cube     = data
        self.orig_cube_hdr = header.copy()
        self.cube          = data.copy()
        self.cube_hdr      = header.copy()
        self.kernel        = nan
        self.zsize         = data.shape[0]
        self.xsize         = data.shape[1]
        self.ysize         = data.shape[2]
        self.vscale        = vel_kms
        self.lam           = lam
        self.camera        = camera
        self.ascale        = scale/1000.

    def rebin(self, new_shape):
        self.cube = zeros((self.zsize, new_shape[0],new_shape[1]))
        self.cube_hdr['CD1_1'] *=  self.orig_cube.shape[1]/new_shape[0]
        self.cube_hdr['CD2_2'] *=  self.orig_cube.shape[2]/new_shape[1]
        self.cube_hdr['CRPIX1'] =  new_shape[0]/2.
        self.cube_hdr['CRPIX2'] =  new_shape[1]/2.
        self.cube_hdr['NAXIS1'] =  new_shape[0]
        self.cube_hdr['NAXIS2'] =  new_shape[1]
        self.cube_hdr['EXTNAME'] =  self.cube_hdr['EXTNAME'].replace('_ORIG', '_OBS')


        for i in arange(self.zsize):
            M, N = self.orig_cube.shape[1:3]
            m, n = new_shape
            if m<M:
                self.cube[i] = self.orig_cube[i].reshape((m,M/m,n,N/n)).mean(3).mean(1)
            else:
                self.cube[i] = np.repeat(np.repeat(self.orig_cube[i], m/M, axis=0), n/N, axis=1)
        self.zsize      = self.cube.shape[0]
        self.xsize      = self.cube.shape[1]
        self.ysize      = self.cube.shape[2]

    def generate_blurred_map(self, kernel_size, band = 'H'):
        self.blrcube    = self.cube.copy()*nan
        self.kernel_size = kernel_size
        self.kernel = Gaussian2DKernel(self.kernel_size)

        print 'Convolving spatially...'
        for i in arange(self.zsize):
            self.blrcube[i] = convolve_fft(self.cube[i], self.kernel)
        #self.blrcube += np.random.normal(0, max(self.blrcube.ravel())/100., shape(self.blrcube))

        #KMOS reaches a point source 5-sigma sensitvity in 8 hr of
        #of (J, H, K) = (22, 21.0, 20.5) AB magnitudes 
        #for R ~ (3380, 3800, 3750)
        #citation: http://www2011.mpe.mpg.de/Highlights/FB2004/exp13_bender.pdf
        if band == 'H': sens, R = 21.0, 3800
        elif band == 'J': sens, R = 22.0, 3380
        elif band == 'K': sens, R = 20.5, 3750
        else: 
            sens = 21.0
            print 'Bad input band, setting sensitivty to H-band value, %s AB mag'%(sens)


        print 'Convolving spectrally...'

        print 'Adding noise...'

        sens_si_fd = (sens*u.ABmag).to(u.Watt/(u.meter*u.meter)/(u.Hz))*astropy.constants.c/(self.lam[0]*u.meter)**2.
        print sens_si_fd

        #The pixel values in our cube are W/m/m^2/Sr (surface brightness). To get to units of 
        #flux density per pixel
        #check to make sure this pixel scale is in proper coordinates
        pix_scale_kpc = self.cube_hdr['CD1_1']*u.kpc
        print pix_scale_kpc, 'per pixel'
        pix_scale_arc = pix_scale_kpc * cosmo.arcsec_per_kpc_proper(1./self.ascale-1)
        print pix_scale_arc, 'per pixel'
        pix_scale_str = (pix_scale_arc**2.).to(u.steradian)
        print pix_scale_str, 'per square pixel'

        #Currently in m, want to get in terms of hz^-1. F_v = (F_lam)*lam^2/c
        #Ths units of this factor are 1/(str*m). Let's now consider an 'aperture' equal to 
        #1 spatial fwhm * 1 spectral fwhm. We want to add noise equal to 1/5th the sensitivity over this aperture.
        #In steradians, the PSF is:
        psf_str = pi*(self.kernel_size**2.*u.arcsec**2.).to(u.steradian)
        print psf_str, 'is the seeing'


        #The spectral lsf fwhm (in pixels) is:
        lsf_pix = (3.e5/R)/(self.vscale[1]-self.vscale[0])
        print lsf_pix
        
        sens_noise = sens_si_fd/psf_str/lsf_pix
        print sens_noise
        self.cube += np.random.normal(0, sens_noise.value, self.cube.shape)

    def generate_intrinsic_kin_map(self):
        self.vel_int    = np.zeros((self.xsize, self.ysize))*np.nan
        self.evel_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.disp_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.edisp_int  = np.zeros((self.xsize, self.ysize))*np.nan
        self.ha_int     = np.zeros((self.xsize, self.ysize))*np.nan

        cb_std = self.orig_cube[0:5,:,:].std()
        
        for i in arange(self.cube.shape[1]):
            for j in arange(self.cube.shape[2]):  
                if self.cube[:,i,j].max() > 10*cb_std: #likely have an Ha line here
                    spec = self.cube[:,i,j]
                    try:
                        c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), self.vscale[nanargmax(spec)], 10, spec[0]])  
                        if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                        (c_a[0] > 0) & (c_a[2] > 5) & (sqrt(v_a[2,2]) < 30) & (sqrt(v_a[1,1]) < 30):
                            self.vel_int[i,j]   = c_a[1]
                            self.evel_int[i,j]  = sqrt(v_a[1,1])
                            self.disp_int[i,j]  = c_a[2]
                            self.edisp_int[i,j] = sqrt(v_a[2,2])
                            self.ha_int[i,j] = c_a[0]*c_a[2]*sqrt(2*pi)                            
                    except:
                        pass

    def generate_observed_kin_map(self):
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
                    if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                    (c_a[0] > 0) & (c_a[2] > 10) & (sqrt(v_a[2,2]) < 30) & (sqrt(v_a[1,1]) < 30):
                        self.vel_obs[i,j]   = c_a[1]
                        self.evel_obs[i,j]  = sqrt(v_a[1,1])
                        self.disp_obs[i,j]  = c_a[2]
                        self.edisp_obs[i,j] = sqrt(v_a[2,2])
                        self.ha_obs[i,j] = c_a[0]*c_a[2]*sqrt(2*pi)                            
                except:
                    pass
        self.vel_obs[-isnan(self.vel_obs)] += np.random.normal(0,5,shape(self.vel_obs[-isnan(self.vel_obs)]))
        self.disp_obs[-isnan(self.disp_obs)] += np.random.normal(0,5,shape(self.vel_obs[-isnan(self.vel_obs)]))

    def get_hdulist(self, master_hdulist):
        colhdr = fits.Header()
        master_hdulist.append(fits.ImageHDU(data = self.orig_cube, header = self.orig_cube_hdr, name = 'cam%i_orig_cube'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.cube, header = self.cube_hdr, name = 'cam%i_obs_cube'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_int, self.edisp_int]), name = 'cam%i_disp_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_obs, self.edisp_obs]), name = 'cam%i_disp_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_int,self.evel_int]), name = 'cam%i_vel_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_obs,self.evel_obs]), name = 'cam%i_vel_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_obs, name = 'cam%i_ha_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_int, name = 'cam%i_ha_int'%self.camera))

        return master_hdulist
