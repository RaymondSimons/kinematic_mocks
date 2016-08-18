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
plt.close('all')
plt.ioff()

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
                spec = self.cube[:,i,j]
                try:
                    c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), 0, 30, spec[0]])           
                    if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                    (c_a[0] > 0) & (sqrt(v_a[2,2]) < 20) & (sqrt(v_a[1,1]) < 20):

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
                spec = self.blrcube[:,i,j]
                try:
                    c_a, v_a = curve_fit(gauss, self.vscale, spec, p0 = [nanmax(spec), 0, 30, spec[0]])           
                    if (isfinite(sqrt(v_a[1,1]))) & (c_a[2] > 0.) & (isfinite(sqrt(v_a[2,2]))) & \
                    (c_a[0] > 0) & (sqrt(v_a[2,2]) < 20) & (sqrt(v_a[1,1]) < 20):

                        self.vel_obs[i,j]   = c_a[1]
                        self.evel_obs[i,j]  = v_a[1,1]
                        self.disp_obs[i,j]  = c_a[2]
                        self.edisp_obs[i,j] = v_a[2,2]

                except:
                    pass

        pass






if __name__ == '__main__':
    mcrx_file = glob.glob('*mcrx.fits')
    if len(mcrx_file) == 0: print 'missing mcrx.fits file'
    else: mcrx_data = fits.open(mcrx_file[0])


    lam = mcrx_data['LAMBDA'].data['lambda']
    c_kms = constants.c.value*1.e-3    #speed of light in km/s
    Ha_m = 6.563e-7 #Halpha in meters
    vel_kms_0 = c_kms*(lam-lam[0])/(lam[0])
    loc_Ha = argmin(abs(lam-Ha_m))
    vel_kms = vel_kms_0 - vel_kms_0[loc_Ha]     
    window = where((vel_kms > -500) & (vel_kms < 500))[0]


    for cam_n in arange(1):
        camera = mcrx_data['CAMERA%i'%(cam_n)]
        kmap = kin_map(camera.data[window], vel_kms[window], lam[window],  cam_n)
        kmap.rebin([40,40])
        kmap.generate_blurred_map(kernel_size = 1)
        kmap.generate_intrinsic_kin_map()
        kmap.generate_observed_kin_map()




        figure(1)
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)

        cmap = matplotlib.cm.jet
        cmap.set_bad('k', 1.0)
        ax1.imshow(kmap.vel_int, cmap = cmap, interpolation = 'nearest')
        ax2.imshow(kmap.disp_int, interpolation = 'nearest')
        ax3.imshow(kmap.vel_obs, cmap = cmap, interpolation = 'nearest')
        ax4.imshow(kmap.disp_obs, interpolation = 'nearest')


        savefig('test.png')
        plt.close('all')














