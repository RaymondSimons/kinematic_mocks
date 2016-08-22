'''
Generate IFU maps from SUNRISE output, need mcrx.fits
created: 08/18/2016, Raymond Simons
'''
import astropy
from astropy.io import fits
import glob
from numpy import ma
import sys
from astropy import constants
import astropy
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel, convolve_fft, convolve
import scipy
from scipy.optimize import curve_fit
import pickle
import pyfits
import os
import numpy as np
from numpy import *

def gauss(x, *p):
    A, mu, sigma, A0= p
    if A < 0:
        return nan*zeros(len(x))
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))+A0


class kin_map():
    def __init__(self, data, vel_kms, lam,  camera):
        self.orig_cube  = data
        self.cube       = data.copy()
        self.kernel     = nan
        self.zsize      = data.shape[0]
        self.xsize      = data.shape[1]
        self.ysize      = data.shape[2]
        self.vscale     = vel_kms
        self.lam        = lam
        self.camera     = camera

    def rebin(self, new_shape):
        self.cube = zeros((self.zsize, new_shape[0],new_shape[1]))

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


    def generate_blurred_map(self, kernel_size):
        self.blrcube    = self.cube.copy()*nan
        self.kernel = Gaussian2DKernel(kernel_size)

        print 'Convolving map...'
        for i in arange(self.zsize):
            self.blrcube[i] = convolve_fft(self.cube[i], self.kernel)

    def generate_intrinsic_kin_map(self):
        self.vel_int    = np.zeros((self.xsize, self.ysize))*np.nan
        self.evel_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.disp_int   = np.zeros((self.xsize, self.ysize))*np.nan
        self.edisp_int  = np.zeros((self.xsize, self.ysize))*np.nan

        for i in arange(self.cube.shape[1]):
            for j in arange(self.cube.shape[2]):  
                krnl = Gaussian1DKernel(1)
                spec = convolve_fft(self.cube[:,i,j], krnl)
                try:
                    c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), 0, 60, spec[0]])  

                    if (sqrt(v_a[1,1]) < 200):
                        self.vel_int[i,j]   = c_a[1]
                        self.evel_int[i,j]  = v_a[1,1]
                        self.disp_int[i,j]  = c_a[2]
                        self.edisp_int[i,j] = v_a[2,2]

                except:
                    pass


    def generate_observed_kin_map(self):
        self.vel_obs    = np.zeros((self.xsize, self.ysize))*np.nan
        self.evel_obs   = np.zeros((self.xsize, self.ysize))*np.nan
        self.disp_obs   = np.zeros((self.xsize, self.ysize))*np.nan
        self.edisp_obs  = np.zeros((self.xsize, self.ysize))*np.nan
        for i in arange(self.blrcube.shape[1]):
            for j in arange(self.blrcube.shape[2]):                
                krnl = Gaussian1DKernel(1)
                spec = convolve_fft(self.blrcube[:,i,j], krnl)

                try:
                    c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), 0, 30, spec[0]])           
                    if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                    (c_a[0] > 0) & (sqrt(v_a[2,2]) < 200) & (sqrt(v_a[1,1]) < 200):

                        self.vel_obs[i,j]   = c_a[1]
                        self.evel_obs[i,j]  = v_a[1,1]
                        self.disp_obs[i,j]  = c_a[2]
                        self.edisp_obs[i,j] = v_a[2,2]

                except:
                    pass

        pass

    def save(self, filename):
        file = open(filename, 'w+')
        file.write(pickle.dumps(self.__dict__))
        file.close()





if __name__ == '__main__':
    kmap_dir = '/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2'


    basepath = os.path.dirname(os.getcwd().replace('/ifu','').replace('_sunrise',''))
    simname = basepath[len(basepath)-6::]
    snapname = os.path.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
    kmap_name = '/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname+'/'+snapname+'/'+snapname
    if not os.path.lexists('/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname):
        os.mkdir('/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname)

    if not os.path.lexists('/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname+'/'+snapname):
        os.mkdir('/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname+'/'+snapname)


    print 'Generating kinematic maps in /nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname+'/'+snapname

    mcrx_file = glob.glob('mcrx.fits')
    kmap_dir = '/nobackupp2/rcsimons/sunrise_testing/kmaps'
    if len(mcrx_file) == 0: print 'missing mcrx.fits file'
    else: mcrx_data = fits.open(mcrx_file[0])


    lam = mcrx_data['LAMBDA'].data['lambda']
    c_kms = constants.c.value*1.e-3    #speed of light in km/s
    Ha_m = 6.563e-7 #Halpha in meters
    vel_kms_0 = c_kms*(lam-lam[0])/(lam[0])
    loc_Ha = argmin(abs(lam-Ha_m))
    vel_kms = vel_kms_0 - vel_kms_0[loc_Ha]     
    window = where((vel_kms > -500) & (vel_kms < 500))[0]

    ncams = mcrx_data['MCRX'].header['N_CAMERA']

    for cam_n in arange(ncams):
        print '\t\t Running on camera %i'%(cam_n)
        camera = mcrx_data['CAMERA%i-NONSCATTER'%(cam_n)]        
        kmap = kin_map(camera.data[window], vel_kms[window], lam[window],  cam_n)
        kmap.rebin([20,20])
        kmap.generate_blurred_map(kernel_size = 0.5)
        kmap.generate_intrinsic_kin_map()
        kmap.generate_observed_kin_map()
        print '\t\t\t Saving to '+kmap_name + '_%i.kmap'%(cam_n)
        kmap.save(kmap_name + '_%i.kmap'%(cam_n))
        

        if False:
            pyfits.writeto('./data/orig_cube.fits', kmap.orig_cube, clobber = True)
            pyfits.writeto('./data/cube.fits', kmap.cube, clobber = True)
            pyfits.writeto('./data/blurred_cube.fits', kmap.blrcube, clobber = True)

        if False:
            fig = figure(1)
            ax1 = fig.add_subplot(221)
            ax2 = fig.add_subplot(222)
            ax3 = fig.add_subplot(223)
            ax4 = fig.add_subplot(224)

            cmap = matplotlib.cm.jet
            cmap.set_bad('k', 1.0)
            ax1.imshow(kmap.vel_int, cmap = cmap, interpolation = 'nearest', vmin = -150, vmax = 150.)
            ax3.imshow(kmap.vel_obs, cmap = cmap, interpolation = 'nearest', vmin = -150, vmax = 150.)

            cmap = matplotlib.cm.viridis
            cmap.set_bad('k', 1.0)
            ax2.imshow(kmap.disp_int, cmap = cmap, interpolation = 'nearest', vmin = 10., vmax = 60.)
            ax4.imshow(kmap.disp_obs, cmap = cmap, interpolation = 'nearest', vmin = 10., vmax = 60.)


            savefig('./figures/test%i.png'%(cam_n))
            plt.close('all')


        if False:
            fig = figure(1)
            ax1 = fig.add_subplot(111)
            cmap = matplotlib.cm.jet
            cmap.set_bad('k', 1.0)
            ax1.plot(kmap.vscale, kmap.orig_cube[:,150,150])
            savefig('./figures/plot%i.png'%(cam_n))
            plt.close('all')













