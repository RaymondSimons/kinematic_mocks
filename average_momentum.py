import numpy as np
import pyfits
from astropy.io import fits
import yt
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
import astropy
from astropy.cosmology import WMAP9 as cosmo
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import scipy
import scipy.ndimage
import scipy.stats as ss
import scipy as sp
from scipy.ndimage.measurements import label
from scipy.optimize import curve_fit
import glob
from scipy.signal import savgol_filter
import os, sys, argparse
from matplotlib.pyplot import *
from numpy import *
from joblib import Parallel, delayed

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    values = values[-isnan(weights)]
    weights = weights[-isnan(weights)]
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))



def write_fits(filename, nir_mstar_cat):
    print '\tGenerating fits for %s...'%filename
    master_hdulist = []
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Storing average momentum measurements in this FITS file."
    prihdr['simname'] = gal
    prihdr['scale'] = aname

    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

    colhdr = fits.Header()

    master_hdulist.append(fits.ImageHDU(data = nir_mstar_cat                                                            , header = colhdr, name = 'nir_mstar_cat'))
    master_hdulist.append(fits.ImageHDU(data = self.L_disk                                                              , header = colhdr, name = 'nir_net_momentum'))
    master_hdulist.append(fits.ImageHDU(data = self.L_disk_s                                                            , header = colhdr, name = 'nir_net_momentum_s'))
    master_hdulist.append(fits.ImageHDU(data = self.stars_id                                                            , header = colhdr, name = 'stars_id'))
    master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_metallicity1 , self.stars_metallicity2))            , header = colhdr, name = 'stars_metallicity'))
    master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_x_cen , self.stars_y_cen , self.stars_z_cen))       , header = colhdr, name = 'stars_xyz_position'))
    master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_vx_cen , self.stars_vy_cen , self.stars_vz_cen))    , header = colhdr, name = 'stars_xyz_velocity'))
    master_hdulist.append(fits.ImageHDU(data = np.stack((self.rr_stars, self.zz_stars))                                 , header = colhdr, name = 'stars_cylindrical_position'))
    master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_jx_cen, self.stars_jy_cen, self.stars_jz_cen))      , header = colhdr, name = 'stars_xyz_momentum'))
    master_hdulist.append(fits.ImageHDU(data = self.epsilon_stars                                                       , header = colhdr, name = 'stars_epsilon'))
    master_hdulist.append(fits.ImageHDU(data = self.mass_profile                                                        , header = colhdr, name = 'mass_profile'))
    master_hdulist.append(fits.ImageHDU(data = self.star_mass                                                           , header = colhdr, name = 'star_mass'))
    master_hdulist.append(fits.ImageHDU(data = self.star_creation_time                                                  , header = colhdr, name = 'star_creation_time'))
    master_hdulist.append(fits.ImageHDU(data = self.star_age                                                            , header = colhdr, name = 'star_age'))

    print '\tSaving to ' + self.fits_name
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_name, clobber = True)

    return


def measure_average_momentum(gal):
    print gal
    '''
    sn_files = glob.glob('../momentum_measurements/%s/*momentum.fits'%(gal))
    anames = array([sn.split('_')[2] for sn in sn_files])
    anames = anames[0:10]
    zs = zeros(len(anames))*nan
    a_arr = zeros(len(anames))*nan 
    mean_jz = zeros((4, len(anames), 2))*nan
    '''





if __name__ == "__main__":

    #args = parse()
    #if args['gal'] is not None: gal = args['gal']
    filename = '/nobackupp2/rcsimons/average_momentum/%s'%gal

    rad_min = 0
    rad_max = 50.
    gals = ['VELA%.2i'%(i+1) for i in np.arange(1)]


    Parallel(n_jobs = -1, backend = 'threading')(delayed(run_momentum_figure)(gal, aname) for gal in gals)





















