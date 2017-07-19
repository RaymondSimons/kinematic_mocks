import astropy
import pyfits
import glob
from glob import glob
import astrodendro
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve_fft
import photutils 
from photutils import detect_sources
from photutils import *
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from joblib import Parallel, delayed
from astropy.io import fits
from numpy import *
import matplotlib.pyplot as plt
import os, sys, argparse
import random

plt.ioff()


def write_fits(fits_name, mom_data, merger_tag, x_stars_box , y_stars_box , z_stars_box, vx_stars_box , vy_stars_box , vz_stars_box):

    print '\tGenerating fits for %s'%fits_name
    master_hdulist = []
    master_hdulist.append(mom_data['PRIMARY'])
    colhdr = fits.Header()
    master_hdulist.append(mom_data['nir_mstar_cat'])
    master_hdulist.append(mom_data['nir_net_momentum'])
    master_hdulist.append(mom_data['nir_net_momentum_s'])
    master_hdulist.append(mom_data['stars_id'])
    master_hdulist.append(fits.ImageHDU(data = np.stack((x_stars_box , y_stars_box , z_stars_box,)), header = colhdr, name = 'stars_xyz_box_position'))
    master_hdulist.append(fits.ImageHDU(data = np.stack((vx_stars_box , vy_stars_box , vz_stars_box)), header = colhdr, name = 'stars_xyz_box_velocity'))
    master_hdulist.append(mom_data['star_mass'])
    master_hdulist.append(mom_data['star_age'])
    master_hdulist.append(fits.ImageHDU(data = merger_tag, header = colhdr, name = 'star_merger_tag'))

    print '\tSaving to ' + fits_name
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_name, clobber = True)

    return master_hdulist


def make_heatmap(ax, epsilon, zz_gas, min_z, max_z, weights = None, good = None, xlabel = 'z height (kpc)', ylabel = 'j$_z$/j$_{circ}$', bins_n = 200, eps_min = 2, eps_max = 2, segm = None, srt_labels = None, do_plot = True):
    if weights == None:
        weights = np.ones(len(zz_gas))

    if good:        
        epsilon = epsilon[good]
        zz_gas = zz_gas[good]
        weights = weights[good]
    heatmap, xedges, yedges = np.histogram2d(epsilon, zz_gas, bins=[linspace(eps_min,eps_max,bins_n), linspace(min_z,max_z,bins_n)], weights = weights)



    sorted_heatmap = argsort(heatmap.ravel())
    vmn = 10.
    vmx_scale = 0.998
    vmx = heatmap.ravel()[sorted_heatmap[int(vmx_scale*len(sorted_heatmap))]]
    heatmap = np.ma.masked_where((heatmap < 10), heatmap)
    heatmap.data[heatmap.data < 10.] = nan

    #heatmap.data[segm > 1] = 0
    if srt_labels!=None:
        #for lbl in srt_labels[1:len(srt_labels)]:
        #    heatmap.data[segm == lbl] = 0
        heatmap.data[segm!=srt_labels[0]] = 0

    if do_plot:
        ax.imshow(heatmap, interpolation = 'nearest', norm = mpl.colors.LogNorm(vmin = vmn, vmax = vmx), origin = 'lower', cmap = 'viridis')
        kern = Gaussian2DKernel(1.)
        kern.normalize()
        heatmap_conv = convolve_fft(heatmap, kern)
        heatmap_conv = np.ma.masked_where((heatmap_conv < 10), heatmap_conv)
        heatmap_conv.data[heatmap_conv.data < 10.] = nan

        X = arange(heatmap.data.shape[0])
        Y = arange(heatmap.data.shape[1])
        Z = log10(heatmap.data)
        ax.contour(X, Y, Z, 8, colors = 'grey')



        ax.set_yticks([0,bins_n/4,bins_n/2,3*bins_n/4,bins_n-1])
        ax.set_xticks([0,bins_n/2,bins_n-1])
        ax.set_xticklabels([format(yedges[0],'.0f'),format(yedges[bins_n/2],'.0f'),format(yedges[bins_n-1],'.0f')])
        ax.set_yticklabels([''])
        ax.set_yticklabels([format(xedges[0],'.0f'),format(xedges[bins_n/4],'.0f'), format(xedges[bins_n/2],'.0f'),format(xedges[3*bins_n/4.],'.0f'),format(xedges[bins_n-1],'.0f')])
        #ax.set_xticklabels([''])
        ax.set_xlabel(xlabel, fontsize = 15)

        ax.set_ylabel(ylabel, fontsize = 20)
        ax.minorticks_on()

        ax.tick_params(axis="both", which='major', color='black', labelcolor='black',size=5, width=1.5)
        ax.tick_params(axis="both", which='minor', color='black', labelcolor='black',size=3, width=1.5)

        return ax, heatmap
    else:
        return heatmap



def add_at(ax, t, loc=2):
    fp = dict(size=10)
    _at = AnchoredText(t, loc=loc, prop=fp)
    ax.add_artist(_at)
    return _at


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    values = values[-isnan(weights)]
    weights = weights[-isnan(weights)]
    average = numpy.average(values, weights=weights)
    variance = numpy.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))


def find_thresh(mn, mx, npix):
    nlabels = 0.
    segm_labels_prev = 0
    mr_prev2 = -99
    mr_prev = -99
    kern = Gaussian2DKernel(0.2, x_size = 4*10, y_size = 4*10)
    kern.normalize()
    a = zeros(kern.array.shape)
    a[kern.array.shape[1]/2.,kern.array.shape[1]/2.] = 1
    kern_2 = Gaussian1DKernel(8)
    a[:,kern.array.shape[1]/2.] = convolve_fft(a[:,kern.array.shape[1]/2.], kern_2)
    a/=sum(a)
    b = convolve_fft(a, kern)
    b/=sum(b)
    temp_heatmap = convolve_fft(heatmap.data, b)
    temp_heatmap[temp_heatmap <= 0] = nan

    for tt, t in enumerate(linspace(mn, mx, 1000)):

        threshold = t
        segm = detect_sources(log10(temp_heatmap), threshold = threshold, npixels = npix)  
        masses = array([sum(temp_heatmap[segm.array == lbl]) for lbl in arange(1, segm.nlabels+1)])
        srt_masses = masses[argsort(masses)[::-1]]
        if len(masses) > 1:
            mass_ratio = srt_masses[0]/srt_masses[1]
            if mr_prev == -99:
                mr_prev = mass_ratio  
                thresh = threshold       
            if (log10(srt_masses[0]) > 7.5) & (log10(srt_masses[1]) > 7.5) & \
                (mr_prev/mass_ratio > 10) & (mass_ratio < 100) & (nansum(srt_masses) > 0.50*nansum(temp_heatmap)):
                thresh = threshold


            mr_prev = mass_ratio

            if len(masses) > 2:
                mass_ratio2 = srt_masses[0]/srt_masses[2]
                if mr_prev2 == -99:
                    mr_prev2 = mass_ratio2
                    thresh = threshold    
                    
                if (log10(srt_masses[0]) > 7.5) & (log10(srt_masses[1]) > 7.5) & (mr_prev2/mass_ratio2 > 10) & (mass_ratio2 < 300) & (nansum(srt_masses) > 0.50*nansum(temp_heatmap)):
                    thresh = threshold
                    
                mr_prev2 = mass_ratio2
        segm_labels_prev = segm.nlabels
    return thresh, temp_heatmap


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

    #parser.add_argument('-s', '--snap_base', default='10MpcBox_csf512_',
    #                    help='Base of the snapshots file names.') 

    #parser.add_argument('-d', '--distance', default=100000, type=float,
    #                    help='Distance between cameras and the center of the galaxy (in [kpc]).')

    #parser.add_argument('--no_export',action='store_true',
    #                    help='Do not export data to fits for Sunrise.') 

    args = vars(parser.parse_args())
    return args




#gals = ['VELA01', 'VELA06', 'VELA07', 'VELA11']

#gals = ['VELA20','VELA21'] 

#gals = ['VELA24', 'VELA27', 'VELA28']

#gals = ['VELA29', 'VELA33', 'VELA34']
#gals = ['VELA34']

#gals = ['VELA01', 'VELA06', 'VELA07', 'VELA11', 'VELA15', 'VELA17', 'VELA20', \
#       'VELA21', 'VELA24', 'VELA27', 'VELA28', 'VELA29', 'VELA33', 'VELA34']

#gals = ['VELA21']
#gals = ['VELA01', 'VELA07', 'VELA11']
#gals = ['VELA01']



#scales = arange(200, 550, 10)
#scales = arange(390, 550, 300)


def run_measure_merger(gal, scale, make_cat = True, do_plot = True):    
    eps_min = -2
    eps_max = 2
    rr_min = 0.
    rr_max = 70
    zz_min = -10
    zz_max = 10
    bins_n = 200

    rec_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/recenter_%s.cat'%gal, skiprows = 12)
    if make_cat: 
        m_cat = open('/nobackupp2/rcsimons/mergers/catalogs/individual/%s_%i.cat'%(gal,scale), 'w+')


    for s, scale in enumerate(scales):
        print gal, '\t', scale
        rec_c = rec_cat[(1000.*rec_cat[:,0]).astype('int') == scale]
        if len(rec_c) > 0:
            rec_c = rec_c[0]
            max_nmergers = 10
            masses_arr = zeros(max_nmergers)*nan
            radii_arr = zeros(max_nmergers)*nan
            jz_arr = zeros(max_nmergers)*nan
            radii_std_arr = zeros(max_nmergers)*nan
            jz_std_arr = zeros(max_nmergers)*nan
            mn_box_pos = zeros((max_nmergers,3))*nan
            mn_box_vel = zeros((max_nmergers,3))*nan


            young_mn = nan
            random.seed(1)
            mom_fl = glob('/nobackupp2/rcsimons/momentum_measurements/%s/*%s*momentum.fits'%(gal, scale))
            rec_fl = glob('/nobackupp2/rcsimons/recenter/%s_%s.fits'%(gal, scale))

            if len(mom_fl) > 0:
                mom_data = fits.open(mom_fl[0])
                rec_data = fits.open(rec_fl[0])
                epsilon_stars = mom_data['STARS_EPSILON'].data
                rr_stars = mom_data['STARS_CYLINDRICAL_POSITION'].data[0]
                zz_stars = mom_data['STARS_CYLINDRICAL_POSITION'].data[1]
                r_stars  = sqrt(sum(mom_data['STARS_XYZ_POSITION'].data**2., axis = 0))

                epsilon_stars_digitized = np.digitize(epsilon_stars, bins = linspace(eps_min, eps_max, bins_n))
                r_stars_digitized = np.digitize(r_stars, bins = linspace(rr_min, rr_max, bins_n))
                empt_arr = np.empty((bins_n-1,bins_n-1), dtype = object)
                for i in arange(bins_n-1):
                    good_r_stars = where(r_stars_digitized == i)[0]
                    r_stars_digitized_new = r_stars_digitized[good_r_stars]
                    epsilon_stars_digitized_new = epsilon_stars_digitized[good_r_stars]                
                    for j in arange(bins_n-1):
                        good_eps_stars = good_r_stars[where(epsilon_stars_digitized_new == j)[0]]
                        empt_arr[i,j] = good_eps_stars


                x_stars_box = rec_data['STARS_XYZ_POSITION_BOX'].data[0]
                y_stars_box = rec_data['STARS_XYZ_POSITION_BOX'].data[1]
                z_stars_box = rec_data['STARS_XYZ_POSITION_BOX'].data[2]

                vx_stars_box = rec_data['STARS_XYZ_VELOCITY_BOX'].data[0]
                vy_stars_box = rec_data['STARS_XYZ_VELOCITY_BOX'].data[1]
                vz_stars_box = rec_data['STARS_XYZ_VELOCITY_BOX'].data[2]

                star_age = mom_data['STAR_AGE'].data
                star_mass= mom_data['STAR_MASS'].data

                if do_plot:
                    plt.close('all')
                    fig  = figure(1, figsize = (25, 5))
                    clf()


                    ax1 = fig.add_subplot(151)
                    ax2 = fig.add_subplot(152)
                    ax3 = fig.add_subplot(153)
                    ax4 = fig.add_subplot(154)
                    ax5 = fig.add_subplot(155)



                    ax1.set_ylabel(r'$\frac{j_z}{j_{circ}}$', fontsize = 30, rotation = 0, labelpad = 20)
                    ax5.set_ylabel(r'$\frac{j_z}{j_{circ}}$', fontsize = 30, rotation = 0, labelpad = 20)

                    rand_arr = np.random.randint(0, len(r_stars), size = 40000)
                    ax1.scatter(r_stars[rand_arr], epsilon_stars[rand_arr], marker = 'o',  s = star_mass[rand_arr]*1.e-3)
                    ax1.set_xlim(rr_min,  rr_max)
                    ax1.set_ylim(eps_min, eps_max)
                    ax1.minorticks_on()
                    ax1.tick_params(axis="both", which='major', color='black', labelcolor='black',size=5, width=1.5)
                    ax1.tick_params(axis="both", which='minor', color='black', labelcolor='black',size=3, width=1.5)


                    ax2, heatmap   = make_heatmap(ax2, epsilon_stars, r_stars, min_z = rr_min, max_z = rr_max, weights = star_mass, 
                                        good = None, xlabel = '', ylabel = '', bins_n = bins_n, eps_min = eps_min, eps_max = eps_max)
                    add_at(ax2, "stars", loc=1)
                else:
                    heatmap   = make_heatmap(None, epsilon_stars, r_stars, min_z = rr_min, max_z = rr_max, weights = star_mass, 
                                        good = None, xlabel = '', ylabel = '', bins_n = bins_n, 
                                        eps_min = eps_min, eps_max = eps_max, do_plot = do_plot)


                npix = 20
                #find_thresh
                mn = 4
                mx = 8
                thresh, temp_heatmap = find_thresh(mn, mx, npix)




                segm = detect_sources(log10(temp_heatmap), threshold = thresh, npixels = npix)     
                m = segm.array
                masked_m = np.ma.masked_where(m == 0, m)
                masses = array([sum(temp_heatmap[segm.array == lbl]) for lbl in arange(1, segm.nlabels+1)])
                st = argsort(masses)[::-1]

                srt_masses = masses[st]
                if sum(srt_masses)/nansum(heatmap.data) < 0.6:
                    mn = 4
                    mx = 6.5
                    thresh, temp_heatmap = find_thresh(mn, mx, npix)
                    segm = detect_sources(log10(temp_heatmap), threshold = thresh, npixels = npix)     
                    m = segm.array
                    masked_m = np.ma.masked_where(m == 0, m)



                if do_plot:
                    pl = ax3.imshow(masked_m, cmap = 'Set1', origin = 'lower', interpolation = 'nearest', vmin = 0., vmax = 8)
                    ax3.set_xticklabels(ax2.get_xticklabels())
                    ax3.set_yticklabels(ax2.get_yticklabels())
                    ax3.set_xticks(ax2.get_xticks())
                    ax3.set_yticks(ax2.get_yticks())
                    ax1.set_xticks([0,35, 70])
                    ax1.set_yticks([-2, -1, 0, 1, 2])

                    ax3.minorticks_on()
                    ax3.tick_params(axis="both", which='major', color='black', labelcolor='black',size=5, width=1.5)
                    ax3.tick_params(axis="both", which='minor', color='black', labelcolor='black',size=3, width=1.5)





                radii = array([weighted_avg_and_std(values = where(segm.array == lbl)[1], weights = temp_heatmap[segm.array == lbl])[0] for lbl in arange(1, segm.nlabels+1)])
                jz = array([weighted_avg_and_std(values = where(segm.array == lbl)[0], weights = temp_heatmap[segm.array == lbl])[0] for lbl in arange(1, segm.nlabels+1)])

                radii_std = array([weighted_avg_and_std(values = where(segm.array == lbl)[1], weights = temp_heatmap[segm.array == lbl])[1] for lbl in arange(1, segm.nlabels+1)])
                jz_std = array([weighted_avg_and_std(values = where(segm.array == lbl)[0], weights = temp_heatmap[segm.array == lbl])[1] for lbl in arange(1, segm.nlabels+1)])


                masses = array([sum(temp_heatmap[segm.array == lbl]) for lbl in arange(1, segm.nlabels+1)])
                st = argsort(masses)[::-1]
                srt_masses = masses[st]
                srt_radii = radii[st]
                srt_radii_std = radii_std[st]
                srt_labels = segm.labels[st]
                srt_jz     = jz[st]
                srt_jz_std  = jz_std[st]

                contours =  segm.outline_segments()
                masked_contours = np.ma.masked_where(contours == 0, contours)

                #plot the correct stars
                merger_tag = np.empty(len(r_stars))

                for i in arange(200-1):
                    for j in arange(200-1):
                        for lll in srt_labels:
                            if masked_m[i,j] == lll:
                                id_list = empt_arr[j,i] #somehow this is swapped, very confused
                                if (id_list != None) & (len(id_list) > 0):
                                    merger_tag[id_list] = lll 
                                    rand_arr = np.random.randint(0, len(id_list), size = min(len(id_list), 1))
                                    id_list = id_list[rand_arr]
                                    ax5.plot(r_stars[id_list], epsilon_stars[id_list], 'k.')

                fits_name = '/nobackupp2/rcsimons/mergers/fits/'+gal+'_a0.'+str(scale)+'_starsmergers.fits'
                master_hdulist = write_fits(fits_name, mom_data, merger_tag, x_stars_box , y_stars_box , z_stars_box, vx_stars_box , vy_stars_box , vz_stars_box)


                mn_box_pos[0:len(masses),0] = array([weighted_avg_and_std(values = x_stars_box[merger_tag == lbl], weights = star_mass[merger_tag == lbl])[0] for lbl in arange(1, segm.nlabels+1)])
                mn_box_pos[0:len(masses),1] = array([weighted_avg_and_std(values = y_stars_box[merger_tag == lbl], weights = star_mass[merger_tag == lbl])[0] for lbl in arange(1, segm.nlabels+1)])
                mn_box_pos[0:len(masses),2] = array([weighted_avg_and_std(values = z_stars_box[merger_tag == lbl], weights = star_mass[merger_tag == lbl])[0] for lbl in arange(1, segm.nlabels+1)])

                mn_box_vel[0:len(masses),0] = array([weighted_avg_and_std(values = vx_stars_box[merger_tag == lbl], weights = star_mass[merger_tag == lbl])[0] for lbl in arange(1, segm.nlabels+1)])
                mn_box_vel[0:len(masses),1] = array([weighted_avg_and_std(values = vy_stars_box[merger_tag == lbl], weights = star_mass[merger_tag == lbl])[0] for lbl in arange(1, segm.nlabels+1)])
                mn_box_vel[0:len(masses),2] = array([weighted_avg_and_std(values = vz_stars_box[merger_tag == lbl], weights = star_mass[merger_tag == lbl])[0] for lbl in arange(1, segm.nlabels+1)])



                ax5.set_xlim(rr_min,  rr_max)
                ax5.set_ylim(eps_min, eps_max)
                ax5.minorticks_on()
                ax5.tick_params(axis="both", which='major', color='black', labelcolor='black',size=5, width=1.5)
                ax5.tick_params(axis="both", which='minor', color='black', labelcolor='black',size=3, width=1.5)
                ax5.set_xticks([0,35, 70])
                ax5.set_yticks([-2, -1, 0, 1, 2])






                if do_plot:
                    #ax2.imshow(masked_contours, cmap = 'Set1', origin = 'lower', vmin = 0., vmax = 8)
                    ax3.annotate(r"%2s%5s%2s%.1f"%('M$_{sum}$','/M$_{tot}$','=',sum(srt_masses)/nansum(heatmap.data)), (107, 55), color = 'black', fontweight = 'bold')

                    if len(masses) > 1:
                        mass_ratio = srt_masses[0]/srt_masses[1]
                        ax3.annotate("%4s%6s%5s"%('m1','',''), (110, 40),  color = cm.Set1(srt_labels[0]/8.), fontweight = 'bold')
                        ax3.annotate("%4s%6s%5s"%('','/m2',''), (110, 40), color = cm.Set1(srt_labels[1]/8.), fontweight = 'bold')
                        ax3.annotate("%4s%6s%5s%.1f"%('','','=',mass_ratio), (110, 40), color = 'black', fontweight = 'bold')
                        ax3.errorbar(srt_radii[0], srt_jz[0], xerr = srt_radii_std[0], yerr = srt_jz_std[0],  fmt = 'o', color = 'black')
                        ax3.errorbar(srt_radii[1], srt_jz[1],  xerr = srt_radii_std[1], yerr = srt_jz_std[1],  fmt = 'o', color = 'black')

                        if len(masses) > 2:
                            mass_ratio = srt_masses[0]/srt_masses[2]
                            ax3.annotate("%4s%6s%5s"%('m1','',''), (110, 25),  color = cm.Set1(srt_labels[0]/8.), fontweight = 'bold')
                            ax3.annotate("%4s%6s%5s"%('','/m3',''), (110, 25), color = cm.Set1(srt_labels[2]/8.), fontweight = 'bold')
                            ax3.annotate("%4s%6s%5s%.1f"%('','','=',mass_ratio), (110, 25), color = 'black', fontweight = 'bold')
                            ax3.errorbar(srt_radii[2], srt_jz[2], xerr = srt_radii_std[2], yerr = srt_jz_std[2],  fmt = 'o', color = 'black')

                            if len(masses) > 3:
                                mass_ratio = srt_masses[0]/srt_masses[3]
                                ax3.annotate("%4s%6s%5s"%('m1','',''), (110, 10),  color = cm.Set1(srt_labels[0]/8.), fontweight = 'bold')
                                ax3.annotate("%4s%6s%5s"%('','/m4',''), (110, 10), color = cm.Set1(srt_labels[3]/8.), fontweight = 'bold')
                                ax3.annotate("%4s%6s%5s%.1f"%('','','=',mass_ratio), (110, 10), color = 'black', fontweight = 'bold')
                                ax3.errorbar(srt_radii[3], srt_jz[3], xerr = srt_radii_std[3], yerr = srt_jz_std[3],  fmt = 'o', color = 'black')

                masses_arr[0:len(masses)] = srt_masses
                radii_arr[0:len(masses)] = srt_radii*(rr_max - rr_min)/temp_heatmap.shape[1] +rr_min 
                jz_arr[0:len(masses)] = srt_jz*(eps_max - eps_min)/temp_heatmap.shape[1] +eps_min 
                radii_std_arr[0:len(masses)] = srt_radii_std*(rr_max - rr_min)/temp_heatmap.shape[1]
                jz_std_arr[0:len(masses)] = srt_jz_std*(eps_max - eps_min)/temp_heatmap.shape[1]

                
                #m = segm.array
                #m_new = convolve_fft(m, kern).astype('int')
                #ax4 = fig.add_subplot(144)
                #masked_mmew = np.ma.masked_where(m_new == 0, m_new)
                #ax4.imshow(masked_mmew, cmap = 'Set1', origin = 'lower', interpolation = 'nearest')
                #ax4.set_xticklabels(ax2.get_xticklabels())
                #ax4.set_yticklabels(ax2.get_yticklabels())
                #ax4.set_xticks(ax2.get_xticks())
                #ax4.set_yticks(ax2.get_yticks())
                #ax4.minorticks_on()
                #ax4.tick_params(axis="both", which='major', color='black', labelcolor='black',size=5, width=1.5)
                #ax4.tick_params(axis="both", which='minor', color='black', labelcolor='black',size=3, width=1.5)


                if do_plot:
                    ax4, heatmap_young   = make_heatmap(ax4, epsilon_stars, r_stars, min_z = rr_min, max_z = rr_max, weights = star_mass, 
                                        good = where(star_age < 20), xlabel = '', ylabel = '', bins_n = bins_n, 
                                        eps_min = eps_min, eps_max = eps_max, segm = segm, srt_labels = srt_labels)
                    ax4.annotate("young stars (<20 Myr)\nof m1", (80, 170), color = 'blue', fontweight = 'bold')

                else:
                    heatmap_young   = make_heatmap(None, epsilon_stars, r_stars, min_z = rr_min, max_z = rr_max, weights = star_mass, 
                                        good = where(star_age < 20), xlabel = '', ylabel = '', bins_n = bins_n, 
                                        eps_min = eps_min, eps_max = eps_max, segm = segm, srt_labels = srt_labels, do_plot = do_plot)







                #for lbl in srt_labels[1:len(srt_labels)]:
                #    heatmap_young[segm.array == lbl] = 0

                #sm = nansum(heatmap_young.data, axis = 1)
                #x = (arange(len(sm))-len(sm)/2.)*(eps_max-eps_min)/(1.*len(sm))
                #young_mn, young_std = weighted_avg_and_std(values = x, weights = sm)

                young_radii, young_radii_std = weighted_avg_and_std(values = where(heatmap_young!= 0)[1], weights = heatmap_young[heatmap_young!= 0])
                young_jz,    young_jz_std    = weighted_avg_and_std(values = where(heatmap_young!= 0)[0], weights = heatmap_young[heatmap_young!= 0])
                if do_plot:
                    ax4.errorbar(young_radii, young_jz, xerr = young_radii_std, yerr = young_jz_std,  fmt = 'o', color = 'black')





                young_rdi_mn     = young_radii*(rr_max - rr_min)/temp_heatmap.shape[1] + rr_min 
                young_rdi_std    = young_radii_std*(rr_max - rr_min)/temp_heatmap.shape[1]

                young_jz_mn     = young_jz*(eps_max - eps_min)/temp_heatmap.shape[1] + eps_min 
                young_jz_std    = young_jz_std*(eps_max - eps_min)/temp_heatmap.shape[1]




                if do_plot:
                    ax1.set_xlabel(r'radius (kpc)', fontsize = 18, rotation = 0, labelpad = 15)
                    ax2.set_xlabel(r'radius (kpc)', fontsize = 18, rotation = 0, labelpad = 15)
                    ax3.set_xlabel(r'radius (kpc)', fontsize = 18, rotation = 0, labelpad = 15)
                    ax4.set_xlabel(r'radius (kpc)', fontsize = 18, rotation = 0, labelpad = 15)
                    ax5.set_xlabel(r'radius (kpc)', fontsize = 18, rotation = 0, labelpad = 15)

                    fig.tight_layout()
                    savefig('/nobackupp2/rcsimons/mergers/figures/merger_maps/%s_%s.png'%(gal, scale), dpi = 300)
                    plt.close('all')
                if make_cat: 
                    #write young
                    m_cat.write('%.3i\t\t'%scale)
                    m_cat.write('%.2f\t'%young_jz_mn)
                    m_cat.write('%.2f\t'%young_jz_std)
                    m_cat.write('%.2f\t'%young_rdi_mn)
                    m_cat.write('%.2f\t'%young_rdi_std)

                    #write all
                    for m, mass in enumerate(masses_arr):
                        if -isnan(mass):
                            m_cat.write('%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t'%(mass/(1.e10), radii_arr[m], radii_std_arr[m], jz_arr[m], jz_std_arr[m], 
                                                                         mn_box_pos[m,0], mn_box_pos[m,1], mn_box_pos[m,2], 
                                                                         mn_box_vel[m,0], mn_box_vel[m,1], mn_box_vel[m,2]))
                            pass
                        else:
                            m_cat.write('%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t'%(mass,  radii_arr[m], radii_std_arr[m],jz_arr[m], jz_std_arr[m],
                                                                         mn_box_pos[m,0], mn_box_pos[m,1], mn_box_pos[m,2], 
                                                                         mn_box_vel[m,0], mn_box_vel[m,1], mn_box_vel[m,2]))
                            pass

                if make_cat: m_cat.write('\n')

    if make_cat: m_cat.close()


if __name__ == "__main__":

    args = parse()
    import yt

    if args['gal'] is not None: gal = args['gal']
    else: print 'no galaxy entered'        
    print "Generating Sunrise Input for: ", gal
    scales = arange(200, 550, 10)
    scales = arange(350, 550, 500)

    Parallel(n_jobs = 1, backend = 'threading')(delayed(run_measure_merger)(gal, scale) for scale in scales)

    m_cat =  open('/nobackupp2/rcsimons/mergers/catalogs/%s.cat'%gal, 'w+')
    cat_hdrs = ['scale',
                'mean jz/jcirc of young stars in central galaxy-- galaxy coordinates',
                'std jz/jcirc of young stars in central galaxy-- galaxy coordinates',
                'mean radial location of young stars in central galaxy (kpc)-- galaxy coordinates',
                'std radial location of young stars in central galaxy (kpc)-- galaxy coordinates',
                'central stellar mass (1.e10 Msun)',
                'central mean radial location (kpc)-- galaxy coordinates',
                'central std radial location  (kpc)-- galaxy coordinates',
                'central mean jz/jcirc-- galaxy coordinates',
                'central std jz/jcirc-- galaxy coordinates',
                'central mean x-position (kpc)-- simulation coordinates',
                'central mean y-position (kpc)-- simulation coordinates',
                'central mean z-position (kpc)-- simulation coordinates',
                'central mean x-velocity (km/s)-- simulation coordinates',
                'central mean y-velocity (km/s)-- simulation coordinates',
                'central mean z-velocity (km/s)-- simulation coordinates',
                'merger 1 stellar mass (1.e10 Msun)',
                'merger 1 mean radial location (kpc)-- galaxy coordinates',
                'merger 1 std radial location  (kpc)-- galaxy coordinates',
                'merger 1 mean jz/jcirc-- galaxy coordinates',
                'merger 1 std jz/jcirc-- galaxy coordinates',
                'merger 1 mean x-position (kpc)-- simulation coordinates',
                'merger 1 mean y-position (kpc)-- simulation coordinates',
                'merger 1 mean z-position (kpc)-- simulation coordinates',
                'merger 1 mean x-velocity (km/s)-- simulation coordinates',
                'merger 1 mean y-velocity (km/s)-- simulation coordinates',
                'merger 1 mean z-velocity (km/s)-- simulation coordinates',
                'merger 1 stellar mass (1.e10 Msun)',
                'merger 2 mean radial location (kpc)-- galaxy coordinates',
                'merger 2 std radial location  (kpc)-- galaxy coordinates',
                'merger 2 mean jz/jcirc-- galaxy coordinates',
                'merger 2 std jz/jcirc-- galaxy coordinates',
                'merger 2 mean x-position (kpc)-- simulation coordinates',
                'merger 2 mean y-position (kpc)-- simulation coordinates',
                'merger 2 mean z-position (kpc)-- simulation coordinates',
                'merger 2 mean x-velocity (km/s)-- simulation coordinates',
                'merger 2 mean y-velocity (km/s)-- simulation coordinates',
                'merger 2 mean z-velocity (km/s)-- simulation coordinates',
                'etc.']

    for i in arange(len(cat_hdrs)):
        if i < len(cat_hdrs):
            m_cat.write('#(%i) %s\n'%(i, cat_hdrs[i]))
        else:
            m_cat.write('#(%i:...) %s\n\n\n\n'%(i, cat_hdrs[i]))
    m_cat.write('\n\n\n\n')

    for s, scale in enumerate(scales):
        cat_s = np.loadtxt('', dtype = 'str', delimiter = 'notarealword')
        m_cat.write('%s\n'%cat_s[0])



    m_cat.close()















