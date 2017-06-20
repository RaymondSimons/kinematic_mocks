import numpy as np
from numpy import *
import os, sys, argparse
import glob
import astropy
from astropy.io import fits
import astropy
from astropy.cosmology import Planck15 as cosmo
from joblib import Parallel, delayed
import scipy
from scipy.interpolate import UnivariateSpline



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

    parser.add_argument('snap_files', nargs='?', default=None, help='Snapshot files to be analyzed.')

    args = vars(parser.parse_args())
    return args




def determine_center(snapfile, nir_cat):
    ds = yt.load(snapfile)
    dd = ds.all_data()
    print 'Loading star velocities...'

    try:
        stars_vx = dd['stars', 'particle_velocity_x'].in_units('km/s')
        assert stars_vx.shape > 5
    except AttributeError,AssertionError:
        print "No star particles found, skipping: ", ds._file_amr
        return 0


    stars_id = dd['stars', 'particle_index']

    stars_vy = dd['stars', 'particle_velocity_y'].in_units('km/s')
    stars_vz = dd['stars', 'particle_velocity_z'].in_units('km/s')

    print 'Loading star positions...'
    stars_x = dd['stars', 'particle_position_x'].in_units('kpc')
    stars_y = dd['stars', 'particle_position_y'].in_units('kpc')
    stars_z = dd['stars', 'particle_position_z'].in_units('kpc')


    print 'Recentering...'
    
    id_cen_star      = nir_cat[1].astype('int')
    cen_star_offset  = nir_cat[2:5].astype('float')
    cen_star_voffset = nir_cat[5:8].astype('float')
    L_disk_s = nir_cat[8:11].astype('float')


    #Determine offset
    cen_x  = stars_x[id_cen_star-1]  - ds.arr(cen_star_offset[0], 'kpc')
    cen_y  = stars_y[id_cen_star-1]  - ds.arr(cen_star_offset[1], 'kpc')
    cen_z  = stars_z[id_cen_star-1]  - ds.arr(cen_star_offset[2], 'kpc')
    cen_vx = stars_vx[id_cen_star-1] - ds.arr(cen_star_voffset[0], 'km/s')
    cen_vy = stars_vy[id_cen_star-1] - ds.arr(cen_star_voffset[1], 'km/s')
    cen_vz = stars_vz[id_cen_star-1] - ds.arr(cen_star_voffset[2], 'km/s')

    return cen_x.value, cen_y.value, cen_z.value, cen_vx.value, cen_vy.value, cen_vz.value, cen_star_offset[0], cen_star_offset[1], cen_star_offset[2], cen_star_voffset[0], cen_star_voffset[1], cen_star_voffset[2] 






if __name__ == "__main__":

    args = parse()
    import yt

    if args['snap_files'] is not None: snaps = [args['snap_files']]
    else: snaps = np.asarray(glob.glob("*.d"))
        
    print "Generating Sunrise Input for: ", snaps

    abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(abssnap)

    dirname = os.path.dirname(abssnap)
    simname = os.path.basename(dirname) #assumes directory name for simulation name

    nir_cat_name = simname[0:-2]+'_v2_'+simname[-2:]
    nir_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Nir_simplified_disc_cat.txt', skiprows = 1, dtype='str')

    print "Simulation name:  ", simname


    new_snapfiles = []
    anames = []
    for sn in snaps:
        aname = sn.split('_')[-1].rstrip('.d')
        snap_dir = os.path.join(simname+'_'+aname+'_sunrise')
        newf = os.path.join(snap_dir,sn)
        anames.append(aname.strip('a'))
        new_snapfiles.append(newf)

    new_snapfiles = np.asarray(new_snapfiles)

    st = ['cen_x', 'cen_y', 'cen_z', 'cen_vx,' 'ven_vy', 'cen_vz', 'cen_star_offset_x', 'cen_star_offset_y', 'cen_star_offset_z', 'cen_star_voffset_x', 'cen_star_voffset_y', 'cen_star_voffset_z']
    out_cat = open('/nobackupp2/rcsimons/catalogs/recenter_%s.cat'%simname, 'w+')
    for l, s in enumerate(st):
        out_cat.write('(%i)%s\n'%(l, s))

    out_cat.write('\n\n\n')
    for i in arange(len(new_snapfiles)):
        print new_snapfiles[i]
        if float(anames[i]) > 0.2:

            in_nir = where(nir_cat[:,0] == 'a'+anames[i])
            print in_nir
            nir_cat_new = nir_cat[in_nir]
            cen_x, cen_y, cen_z, cen_vx, cen_vy, cen_vz, cen_star_offset_x, cen_star_offset_y, cen_star_offset_z, cen_star_voffset_x, cen_star_voffset_y, cen_star_voffset_z  = determine_center(new_snapfiles[i], nir_cat_new)
            out_cat.write('%5s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n'%(anames[i], cen_x, cen_y, cen_z, cen_vx, cen_vy, cen_vz, cen_star_offset_x, cen_star_offset_y, cen_star_offset_z, cen_star_voffset_x, cen_star_voffset_y, cen_star_voffset_z))



    out_cat.close()

















