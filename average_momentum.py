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



def write_fits(fits_filename, gal, zs, mean_jz, smoothed_mn):
    print '\tGenerating fits for %s...'%fits_filename
    master_hdulist = []
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Storing average momentum measurements in this FITS file."
    prihdr['simname'] = gal

    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

    colhdr = fits.Header()

    master_hdulist.append(fits.ImageHDU(data = zs, header = colhdr, name = 'redshift'))
    master_hdulist.append(fits.ImageHDU(data = mean_jz, header = colhdr, name = 'mean_jz'))
    master_hdulist.append(fits.ImageHDU(data = smoothed_mn , header = colhdr, name = 'smoothed_mn'))

    print '\tSaving to ' + fits_filename
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_filename, clobber = True)

    return


def measure_average_momentum(gal):
    print gal
    fits_filename = '/nobackupp2/rcsimons/average_momentum/%s_avg_momentum.fits'%gal
    sn_files = glob.glob('/nobackupp2/rcsimons/momentum_measurements/%s/*momentum.fits'%(gal))
    anames = array([sn.split('_')[2] for sn in sn_files])
    zs = zeros(len(anames))*nan
    a_arr = zeros(len(anames))*nan 
    mean_jz = zeros((4, len(anames), 2))*nan
    for a, aname in enumerate(anames):
        a_arr[a] = float(aname.strip('a'))
        zs[a] = 1./a_arr[a] - 1.
        print '\t', aname, zs[a]
        data = pyfits.open('/nobackupp2/rcsimons/momentum_measurements/%s/%s_%s_momentum.fits'%(gal, gal, aname))
        epsilon_stars = data['STARS_EPSILON'].data
        rad_stars = sqrt(sum(data['STARS_XYZ_POSITION'].data**2., axis = 0))
        star_age=data['STAR_AGE'].data
        star_mass=data['STAR_MASS'].data

        gas_epsilon = sum(data['GAS_RAD_EPSILON'].data[:,0:100], axis = 1) #only summing up within 50 kpc
        gas_epsilon_edges =  data['GAS_RAD_EPSILON_EDGES'].data[0]
        min_eps = min(gas_epsilon_edges)
        max_eps = max(gas_epsilon_edges)

        mn, st = weighted_avg_and_std(arange(len(gas_epsilon)), gas_epsilon)
        mn = mn*(max_eps - min_eps)/len(gas_epsilon) + min_eps
        st = st*(max_eps - min_eps)/len(gas_epsilon)

        mean_jz[0,a,0], mean_jz[0,a,1] = mn, st


        good_young = where((rad_stars < rad_max) & (star_age < 20.e6) & isfinite(epsilon_stars))[0]
        good_intermediate = where((rad_stars < rad_max) & (star_age > 1.e8) & (star_age < 3.e8) & isfinite(epsilon_stars))[0]
        good_old = where((rad_stars < rad_max) & (star_age > 1.e9) & isfinite(epsilon_stars))[0]

        for g, good in enumerate(array([good_young, good_intermediate, good_old])):
            if len(good) > 0:
                mn, st = weighted_avg_and_std(epsilon_stars[good], star_mass[good])
                mean_jz[g+1,a,0], mean_jz[g+1,a,1] = mn, st

    smoothed_mn = zeros((4, len(anames), 3))

    wl = 5
    po = 0

    sorted_zs = argsort(zs)

    mean_jz = mean_jz[:,sorted_jz,:]
    zs = zs[sorted_jz]



    for i in arange(len(smoothed_mn)):
        print i, mean_jz.shape
        good = where(-isnan(mean_jz[i,:,0]))
        smoothed_mn[i,good,0] = savgol_filter(mean_jz[i,good,0], window_length=wl, polyorder=po)
        smoothed_mn[i,good,1] = savgol_filter(mean_jz[i,good,0]+mean_jz[i,good,1], window_length=wl, polyorder=po)
        smoothed_mn[i,good,2] = savgol_filter(mean_jz[i,good,0]-mean_jz[i,good,1], window_length=wl, polyorder=po)

    write_fits(fits_filename, gal, zs, mean_jz, smoothed_mn)
    return


if __name__ == "__main__":

    #args = parse()
    #if args['gal'] is not None: gal = args['gal']

    rad_min = 0
    rad_max = 50.
    gals = ['VELA%.2i'%(i+1) for i in np.arange(35)]


    Parallel(n_jobs = -1, backend = 'threading')(delayed(measure_average_momentum)(gal) for gal in gals)





















