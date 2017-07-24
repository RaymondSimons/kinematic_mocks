import numpy as np
import pyfits
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
import os, sys, argparse
from matplotlib.pyplot import *
from numpy import *
from joblib import Parallel, delayed



def make_heatmap(ax, epsilon, zz_gas, min_z, max_z, weights = None, good = None, xlabel = 'z height (kpc)', ylabel = 'j$_z$/j$_{circ}$', bins_n = 200, eps_min = 2, eps_max = 2):
    if weights == None:
        weights = np.ones(len(zz_gas))

    if good:        
        epsilon = epsilon[good]
        zz_gas = zz_gas[good]
        weights = weights[good]
    heatmap, xedges, yedges = np.histogram2d(epsilon, zz_gas, bins=[linspace(eps_min,eps_max,bins_n), linspace(min_z,max_z,bins_n)], weights = weights)
    sorted_heatmap = argsort(heatmap.ravel())
    vmn = 0.
    vmx_scale = 0.998
    vmx = heatmap.ravel()[sorted_heatmap[int(vmx_scale*len(sorted_heatmap))]]
    heatmap = np.ma.masked_where((0. == heatmap), heatmap)

    ax.imshow(heatmap, interpolation = 'nearest',vmin = vmn, vmax = vmx, origin = 'lower', cmap = 'viridis')

    ax.set_yticks([0,bins_n/4,bins_n/2,3*bins_n/4,bins_n-1])
    ax.set_xticks([0,bins_n/2,bins_n-1])
    ax.set_xticklabels([format(yedges[0],'.0f'),format(yedges[bins_n/2],'.0f'),format(yedges[bins_n-1],'.0f')])
    ax.set_yticklabels([''])
    ax.set_yticklabels([format(xedges[0],'.0f'),format(xedges[bins_n/4],'.0f'), format(xedges[bins_n/2],'.0f'),format(xedges[3*bins_n/4.],'.0f'),format(xedges[bins_n-1],'.0f')])
    #ax.set_xticklabels([''])
    ax.set_xlabel(xlabel, fontsize = 15)

    ax.set_ylabel(ylabel, fontsize = 30, rotation = 0, labelpad = 24)
    ax.minorticks_on()

    ax.tick_params(axis="both", which='major', color='black', labelcolor='black',labelsize=10, size=5, width=1.5)
    ax.tick_params(axis="both", which='minor', color='black', labelcolor='black',labelsize=10, size=3, width=1.5)

    return ax


def make_cold_gas_heatmap(ax, heatmap, xedges, yedges, xlabel = 'z height (kpc)', ylabel = 'j$_z$/j$_{circ}$'):
    bins_n = len(xedges)
    sorted_heatmap = argsort(heatmap.ravel())
    vmn = 0.
    vmx_scale = 0.998
    vmx = heatmap.ravel()[sorted_heatmap[int(vmx_scale*len(sorted_heatmap))]]
    heatmap = np.ma.masked_where((0. == heatmap), heatmap)

    ax.imshow(heatmap, interpolation = 'nearest',vmin = vmn, vmax = vmx, origin = 'lower', cmap = 'viridis')

    ax.set_yticks([0,bins_n/4,bins_n/2,3*bins_n/4,bins_n-1])
    ax.set_xticks([0,bins_n/2,bins_n-1])
    ax.set_xticklabels([format(yedges[0],'.0f'),format(yedges[bins_n/2],'.0f'),format(yedges[bins_n-1],'.0f')])
    ax.set_yticklabels([''])
    ax.set_yticklabels([format(xedges[0],'.0f'),format(xedges[bins_n/4],'.0f'), format(xedges[bins_n/2],'.0f'),format(xedges[3*bins_n/4.],'.0f'),format(xedges[bins_n-1],'.0f')])
    ax.set_xlabel(xlabel, fontsize = 15)

    ax.set_ylabel(ylabel, fontsize = 30, rotation = 0, labelpad = 24)
    ax.minorticks_on()



    ax.tick_params(axis="both", which='major', color='black', labelcolor='black',labelsize=10, size=5, width=1.5)
    ax.tick_params(axis="both", which='minor', color='black', labelcolor='black',labelsize=10, size=3, width=1.5)

    return ax




def add_at(ax, t, loc=2):
    fp = dict(size=10)
    _at = AnchoredText(t, loc=loc, prop=fp)
    ax.add_artist(_at)
    return _at


#This file will be used to store the profile of the momentum
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

    parser.add_argument('gal', nargs='?', default=None, help='Galaxy to be analyzed')

    args = vars(parser.parse_args())
    return args

    def run_momentum_figure(gal, aname):
        print 'Making plot...'
        eps_min = -2
        eps_max = 2
        rr_min = 0
        rr_max = 30
        zz_min = -10
        zz_max = 10
        rad_min = 0
        rad_max = 100.

        print aname
        plt.ioff()
        plt.close('all')
        fig = figure(figsize=(9,12))

        a = pyfits.open('/nobackupp2/rcsimons/momentum_measurements/%s/%s_%s_momentum.fits'%(gal, gal, aname))

        epsilon_stars = a['STARS_EPSILON'].data
        rr_stars=a['STARS_CYLINDRICAL_POSITION'].data[0]
        zz_stars=a['STARS_CYLINDRICAL_POSITION'].data[1]
        rad_stars = sqrt(sum(a['STARS_XYZ_POSITION'].data**2., axis = 0))

        star_age=a['STAR_AGE'].data
        star_mass=a['STAR_MASS'].data

        ax = fig.add_subplot(431)
        cg_str = "Cold Gas\n"+r"T < 10$^4$ K"
        ax = make_cold_gas_heatmap(ax, a['GAS_ZZ_EPSILON'].data, a['GAS_ZZ_EPSILON_EDGES'].data[0], a['GAS_ZZ_EPSILON_EDGES'].data[1],
                                     xlabel = '', ylabel = r'$\frac{j_z}{j_{circ}}$')
        add_at(ax, cg_str, loc=4)


        ax = fig.add_subplot(432)
        ax = make_cold_gas_heatmap(ax, a['GAS_RR_EPSILON'].data, a['GAS_RR_EPSILON_EDGES'].data[0], a['GAS_RR_EPSILON_EDGES'].data[1],
                                     xlabel = '', ylabel = '')
        add_at(ax, cg_str, loc=4)

        ax.set_title(gal+"\n"+r"$z=%.2f$"%(1./float(aname.strip('a'))-1.), fontweight = 'bold')

        ax = fig.add_subplot(433)
        ax = make_cold_gas_heatmap(ax, a['GAS_RAD_EPSILON'].data, a['GAS_RAD_EPSILON_EDGES'].data[0], a['GAS_RAD_EPSILON_EDGES'].data[1],
                                     xlabel = '', ylabel = '')
        add_at(ax, cg_str, loc=4)





        ax = fig.add_subplot(434)
        good = where((abs(rr_stars) < rr_max) & (star_age < 20.e6))
        ys_str = "Young stars\nage < 20 Myr"
        ax = make_heatmap(ax, epsilon_stars, zz_stars, zz_min, zz_max, weights = star_mass, good = good, xlabel = '', ylabel = r'$\frac{j_z}{j_{circ}}$', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, ys_str, loc=4)
        

        ax = fig.add_subplot(435)
        good = where((abs(zz_stars) < zz_max) & (star_age < 20.e6) & isfinite(epsilon_stars))
        ax = make_heatmap(ax, epsilon_stars, rr_stars, rr_min, rr_max, weights = star_mass, good = good, xlabel = '', ylabel = '', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, ys_str, loc=4)


        ax = fig.add_subplot(436)
        good = where((rad_stars < rad_max) & (star_age < 20.e6) & isfinite(epsilon_stars))
        ax = make_heatmap(ax, epsilon_stars, rad_stars, rad_min, rad_max, weights = star_mass, good = good, xlabel = '', ylabel = '', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, ys_str, loc=4)



        is_str = "Intermediate stars\n100 < age < 300 Myr"
        ax = fig.add_subplot(437)
        good = where((abs(rr_stars) < rr_max) & (star_age > 1.e8) & (star_age < 3.e8))
        ax = make_heatmap(ax, epsilon_stars, zz_stars, zz_min, zz_max, weights = star_mass, good = good, xlabel = '', ylabel = r'$\frac{j_z}{j_{circ}}$', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, is_str, loc=4)
        

        ax = fig.add_subplot(438)
        good = where((abs(zz_stars) < zz_max) & (star_age > 1.e8) & (star_age < 3.e8) & isfinite(epsilon_stars))
        ax = make_heatmap(ax, epsilon_stars, rr_stars, rr_min, rr_max, weights = star_mass, good = good, xlabel = '', ylabel = '', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, is_str, loc=4)

        ax = fig.add_subplot(439)
        good = where((rad_stars < rad_max) & (star_age > 1.e8) & (star_age < 3.e8) & isfinite(epsilon_stars))
        ax = make_heatmap(ax, epsilon_stars, rad_stars, rad_min, rad_max, weights = star_mass, good = good, xlabel = '', ylabel = '', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, is_str, loc=4)




        ax = fig.add_subplot(4,3,10)
        os_str = "Old stars\nage > 1 Gyr"
        good = where((abs(rr_stars) < rr_max) & (star_age > 1.e9))
        ax = make_heatmap(ax, epsilon_stars, zz_stars, zz_min, zz_max, weights = star_mass, good = good, xlabel = 'distance above disk\n(kpc)', ylabel = r'$\frac{j_z}{j_{circ}}$', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, os_str, loc=4)



        ax = fig.add_subplot(4,3,11)
        good = where((abs(zz_stars) < zz_max) & (star_age > 1.e9) & isfinite(epsilon_stars))
        ax = make_heatmap(ax, epsilon_stars, rr_stars, rr_min, rr_max, weights = star_mass, good = good, xlabel = 'distance along disk\n(kpc)', ylabel = '', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)

        add_at(ax, os_str, loc=4)

        ax = fig.add_subplot(4,3,12)
        good = where((rad_stars < rad_max) & (star_age > 1.e9) & isfinite(epsilon_stars))
        ax = make_heatmap(ax, epsilon_stars, rad_stars, rad_min, rad_max, weights = star_mass, good = good, xlabel = 'distance radial\n(kpc)', ylabel = '', 
                     bins_n = 50, eps_min = eps_min, eps_max = eps_max)
        add_at(ax, os_str, loc=4)



        fig.subplots_adjust(hspace = 0.0)

        print 'Saving plot...'

        savefig('/nobackupp2/rcsimons/figures/momentum_figures/%s_%s_momentum_heat.png'%(gal, aname), dpi = 300)
        plt.close('all')
        return


if __name__ == "__main__":

    args = parse()
    import yt
    if args['gal'] is not None: gal = args['gal']

    sn_files = glob.glob('/nobackupp2/rcsimons/momentum_measurements/%s/*momentum.fits'%(gal))

    anames = array([sn.split('_')[2] for sn in sn_files])
    anames = anames

    Parallel(n_jobs = -1, backend = 'threading')(delayed(run_momentum_figure)(gal, aname) for aname in anames)




