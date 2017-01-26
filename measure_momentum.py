import numpy as np
from numpy import *
import os, sys, argparse
import glob
import astropy
from astropy.io import fits
import astropy
from astropy.cosmology import Planck15 as cosmo
from joblib import Parallel, delayed



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

    #parser.add_argument('-s', '--snap_base', default='10MpcBox_csf512_',
    #                    help='Base of the snapshots file names.') 

    #parser.add_argument('-d', '--distance', default=100000, type=float,
    #                    help='Distance between cameras and the center of the galaxy (in [kpc]).')

    #parser.add_argument('-f', '--fov', default=50, type=float,
    #                    help='Field of view of the cameras at the image plane (in [kpc]).')

    #parser.add_argument('--no_export',action='store_true',
    #                    help='Do not export data to fits for Sunrise.') 

    args = vars(parser.parse_args())
    return args





class momentum_obj():
    def __init__(self, simname, aname, snapfile, fits_name):
        self.simname = simname
        self.aname = aname
        self.snapfile = snapfile
        self.fits_name = 


    def write(self):
        print self.simname
        print self.aname
        print self.snapfile


    def load(self):
        ds = yt.load(self.snapfile)
        dd = ds.all_data()

        print 'Loading gas velocity...'
        self.gas_vx = dd['gas', 'velocity_x']
        self.gas_vy = dd['gas', 'velocity_y']
        self.gas_vz = dd['gas', 'velocity_z']

        print 'Loading gas temperature...'
        self.gas_temp = dd['gas', 'temperature']

        print 'Loading gas cell mass...'
        self.gas_mass = dd['gas', 'cell_mass']

        print 'Loading cell potential...'
        self.gas_potential = dd['gas', 'potential']

        print 'Loading gas cell position...'
        self.gas_x = dd['index', 'x'].in_units('kpc')
        self.gas_y = dd['index', 'y'].in_units('kpc')
        self.gas_z = dd['index', 'z'].in_units('kpc')



        print 'Loading star positions...'
        self.star_x = dd['stars', 'particle_position_x'].in_units('kpc')
        self.star_y = dd['stars', 'particle_position_y'].in_units('kpc')
        self.star_z = dd['stars', 'particle_position_z'].in_units('kpc')
        self.star_mass = dd['stars', 'particle_mass'].in_units('Msun')
        self.star_creation_time = dd['stars', 'particle_creation_time'].in_units('yr')
        self.star_age = ds.arr(cosmo.age(ds.current_redshift).value, 'Gyr').in_units('yr') - self.star_creation_time

        print 'Loading star velocities...'
        self.star_vx = dd['stars', 'particle_velocity_x'].in_units('km/s')
        self.star_vy = dd['stars', 'particle_velocity_y'].in_units('km/s')
        self.star_vz = dd['stars', 'particle_velocity_z'].in_units('km/s')



    def write_fits(self):
        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the momentum measurements in this FITS file."
        prihdr['simname'] = self.simname
        prihdr['scale'] = self.aname.strip('a')
        prihdr['snapfile'] = self.snapfile

        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)


        '''
        colhdr = fits.Header()
        master_hdulist.append(fits.ImageHDU(data = self.orig_cube, header = self.orig_cube_hdr, name = 'cam%i_orig_cube'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.cube, header = self.cube_hdr, name = 'cam%i_obs_cube'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_int, self.edisp_int]), name = 'cam%i_disp_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_obs, self.edisp_obs]), name = 'cam%i_disp_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_int,self.evel_int]), name = 'cam%i_vel_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_obs,self.evel_obs]), name = 'cam%i_vel_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_obs, name = 'cam%i_ha_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_int, name = 'cam%i_ha_int'%self.camera))
        '''

        thdulist = fits.HDUList(col_list)
        thdulist.writeto(out_dir+'/'+'kinematics.fits', clobber = True)        


        thdulist = fits.HDUList(master_hdulist)
        print '\t\t\t Saving to ' + self.fits_name
        thdulist.writeto(self.fits_name, clobber = True)



        return master_hdulist







def measure_momentum(snapfile):
    print 'Measuring momentum for '+ snapfile
    aname = (os.path.basename(snapfile)).split('_')[-1].rstrip('.d')
    simname = snapfile.split('_')[0]
    fits_name = simname+'_'+aname+'_momentum.fits'
    mom = momentum_obj(simname, aname, snapfile, fits_name)
    mom.write()
    #mom.load()
    mom.write_fits()





if __name__ == "__main__":

    args = parse()
    import yt
    print args

    if args['snap_files'] is not None:
        snaps = [args['snap_files']]
    else:
        snaps = np.asarray(glob.glob("*.d"))

        
    print "Generating Sunrise Input for: ", snaps

    abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(abssnap)

    dirname = os.path.dirname(abssnap)
    simname = os.path.basename(dirname) #assumes directory name for simulation name
    print "Simulation name:  ", simname

    out_sim_dir = os.path.join('/nobackupp2/rcsimons/momentum_measurements/', simname)
    print os.path.lexists(out_sim_dir)
    if not os.path.lexists(out_sim_dir):
        print 'Creating momentum directory for %s'%out_sim_dir
        os.mkdir(out_sim_dir)                



    particle_headers = []
    particle_data = []
    stars_data = []
    new_snapfiles = []
    for sn in snaps:
        aname = sn.split('_')[-1].rstrip('.d')
        particle_headers.append('PMcrd'+aname+'.DAT')
        particle_data.append('PMcrs0'+aname+'.DAT')
        stars_data.append('stars_'+aname+'.dat')
        snap_dir = os.path.join(simname+'_'+aname+'_sunrise')
        
        #print "Sunrise directory: ", snap_dir
        if not os.path.lexists(snap_dir):
            os.mkdir(snap_dir)        

        newf = os.path.join(snap_dir,sn)
        new_snapfiles.append(newf)
        if not os.path.lexists(newf):
            os.symlink(os.path.abspath(sn),newf)
            os.symlink(os.path.abspath(particle_headers[-1]),os.path.join(snap_dir,particle_headers[-1]))
            os.symlink(os.path.abspath(particle_data[-1]),os.path.join(snap_dir,particle_data[-1]))
            os.symlink(os.path.abspath(stars_data[-1]),os.path.join(snap_dir,stars_data[-1]))


    new_snapfiles = np.asarray(new_snapfiles)


    Parallel(n_jobs = -1)(delayed(measure_momentum)(new_snapfiles[i]) for i in arange(len(new_snapfiles)))



    '''
    # Loop over snapshot to generate cameras and projection plots, 
    # parallelization happens while generating the plots.
    galprops_file = simname+'_galprops.npy'
    galprops = np.load(galprops_file)[()]



    for snapfile in new_snapfiles:

        if os.path.abspath(snapfile) not in galprops['snap_files']: continue
        idx = np.argwhere(galprops['snap_files']==os.path.abspath(snapfile))[0][0]


        ds = yt.load(snapfile)
        dd = ds.all_data()

        print 'Loading gas velocity...'
        gas_vx = dd['gas', 'velocity_x']
        gas_vy = dd['gas', 'velocity_y']
        gas_vz = dd['gas', 'velocity_z']

        print 'Loading gas temperature...'
        gas_temp = dd['gas', 'temperature']

        print 'Loading gas cell mass...'
        gas_mass = dd['gas', 'cell_mass']

        print 'Loading cell potential...'
        gas_potential = dd['gas', 'potential']

        print 'Loading gas cell position...'
        gas_x = dd['index', 'x'].in_units('kpc')
        gas_y = dd['index', 'y'].in_units('kpc')
        gas_z = dd['index', 'z'].in_units('kpc')



        print 'Loading star positions...'
        star_x = dd['stars', 'particle_position_x'].in_units('kpc')
        star_y = dd['stars', 'particle_position_y'].in_units('kpc')
        star_z = dd['stars', 'particle_position_z'].in_units('kpc')
        star_mass = dd['stars', 'particle_mass'].in_units('Msun')
        star_creation_time = dd['stars', 'particle_creation_time'].in_units('yr')
        star_age_all = ds.arr(cosmo.age(ds.current_redshift).value, 'Gyr').in_units('yr') - star_creation_time

        print 'Loading star velocities...'
        star_vx = dd['stars', 'particle_velocity_x'].in_units('km/s')
        star_vy = dd['stars', 'particle_velocity_y'].in_units('km/s')
        star_vz = dd['stars', 'particle_velocity_z'].in_units('km/s')


        stars_L = galprops['stars_L'][idx]
        gas_L 	= galprops['gas_L'][idx]


        try:
            L_sum = stars_L + gas_L
        except TypeError:
            L_sum = gas_L

            
        L = L_sum/np.sqrt(np.sum(L_sum*L_sum))

        print stars_L, gas_L



        col_list = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the momentum properties in this FITS file."
        prihdu = fits.PrimaryHDU(header=prihdr)
        col_list.append(prihdu)
        col_list.append(fits.ImageHDU(data = young_stars_L, name = 'Young_Stars_L'))        
        col_list.append(fits.ImageHDU(data = stars_L, name = 'Stars_L'))
        col_list.append(fits.ImageHDU(data = gas_L, name = 'Gas_L'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Gas_circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Gas_circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Gas_circularity_r'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Young_Stars_circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Young_Stars_circularity_r'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Old_Stars_circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Old_Stars_circularity_r'))

        thdulist = fits.HDUList(col_list)
        thdulist.writeto(out_dir+'/'+'kinematics.fits', clobber = True)
    '''
















