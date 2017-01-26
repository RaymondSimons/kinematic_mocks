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
        self.fits_name = fits_name


    def write(self):
        print self.simname
        print self.aname
        print self.snapfile


    def load(self):
        ds = yt.load(self.snapfile)
        dd = ds.all_data()


        print 'Loading star velocities...'
        self.stars_vx = dd['stars', 'particle_velocity_x'].in_units('km/s')
        self.stars_vy = dd['stars', 'particle_velocity_y'].in_units('km/s')
        self.stars_vz = dd['stars', 'particle_velocity_z'].in_units('km/s')

        print 'Loading star positions...'
        self.stars_x = dd['stars', 'particle_position_x'].in_units('kpc')
        self.stars_y = dd['stars', 'particle_position_y'].in_units('kpc')
        self.stars_z = dd['stars', 'particle_position_z'].in_units('kpc')

        print 'Loading gas velocity...'
        self.gas_vx = dd['gas', 'velocity_x'].in_units('km/s')
        self.gas_vy = dd['gas', 'velocity_y'].in_units('km/s')
        self.gas_vz = dd['gas', 'velocity_z'].in_units('km/s')
        

        print 'Loading gas cell position...'
        self.gas_x = dd['index', 'x'].in_units('kpc')
        self.gas_y = dd['index', 'y'].in_units('kpc')
        self.gas_z = dd['index', 'z'].in_units('kpc')


        '''
        print 'Loading gas temperature...'
        self.gas_temp = dd['gas', 'temperature']

        print 'Loading gas cell mass...'
        self.gas_mass = dd['gas', 'cell_mass']

        print 'Loading cell potential...'
        self.gas_potential = dd['gas', 'potential']




        self.star_mass = dd['stars', 'particle_mass'].in_units('Msun')
        self.star_creation_time = dd['stars', 'particle_creation_time'].in_units('yr')
        self.star_age = ds.arr(cosmo.age(ds.current_redshift).value, 'Gyr').in_units('yr') - self.star_creation_time

        '''

        print 'Finished loading...'


    def calc_momentum(self, nir_cat, nir_disc_cat):
        print 'Calculating momentum...'

        ds = yt.load(self.snapfile, limit_level = 8)
        
        id_cen_star      = nir_cat[1].astype('int')
        cold_cen         = nir_disc_cat[1:4].astype('float')
        cen_star_offset  = nir_cat[2:5].astype('float')
        cen_star_voffset = nir_cat[5:8].astype('float')

        #Determine offset
        cen_x  = self.stars_x[id_cen_star-1]  - ds.arr(cen_star_offset[0], 'kpc')
        cen_y  = self.stars_y[id_cen_star-1]  - ds.arr(cen_star_offset[1], 'kpc')
        cen_z  = self.stars_z[id_cen_star-1]  - ds.arr(cen_star_offset[2], 'kpc')
        cen_vx = self.stars_vx[id_cen_star-1] - ds.arr(cen_star_voffset[0], 'km/s')
        cen_vy = self.stars_vy[id_cen_star-1] - ds.arr(cen_star_voffset[1], 'km/s')
        cen_vz = self.stars_vz[id_cen_star-1] - ds.arr(cen_star_voffset[2], 'km/s')




        #Recenter positions and velocities for stars
        self.stars_x_cen   = self.stars_x  - cen_x
        self.stars_y_cen   = self.stars_y  - cen_y
        self.stars_z_cen   = self.stars_z  - cen_z
        self.stars_pos_mag = sqrt(self.stars_x_cen**2.  + self.stars_y_cen**2.  + self.stars_z_cen**2.)

        self.stars_vx_cen  = self.stars_vx - cen_vx
        self.stars_vy_cen  = self.stars_vy - cen_vy
        self.stars_vz_cen  = self.stars_vz - cen_vz
        self.stars_vel_mag = sqrt(self.stars_vx_cen**2. + self.stars_vy_cen**2. + self.stars_vz_cen**2.)

        #Recenter positions and velocities for gas
        self.gas_x_cen   = self.gas_x  - cen_x
        self.gas_y_cen   = self.gas_y  - cen_y
        self.gas_z_cen   = self.gas_z  - cen_z
        self.gas_pos_mag = sqrt(self.gas_x_cen**2.  + self.gas_y_cen**2.  + self.gas_z_cen**2.)

        self.gas_vx_cen  = self.gas_vx - cen_vx
        self.gas_vy_cen  = self.gas_vy - cen_vy
        self.gas_vz_cen  = self.gas_vz - cen_vz
        self.gas_vel_mag = sqrt(self.gas_vx_cen**2. + self.gas_vy_cen**2. + self.gas_vz_cen**2.)





        #Calculate momentum for stars
        self.stars_jx_cen = self.stars_vz_cen * self.stars_y_cen - self.stars_z_cen * self.stars_vy_cen
        self.stars_jy_cen = self.stars_vx_cen * self.stars_z_cen - self.stars_x_cen * self.stars_vz_cen
        self.stars_jz_cen = self.stars_vy_cen * self.stars_x_cen - self.stars_y_cen * self.stars_vx_cen
        self.stars_j_mag  = sqrt(self.stars_jx_cen**2. + self.stars_jy_cen**2. + self.stars_jz_cen**2.)


        #Calculate momentum for gas
        self.gas_jx_cen = self.gas_vz_cen * self.gas_y_cen - self.gas_z_cen * self.gas_vy_cen
        self.gas_jy_cen = self.gas_vx_cen * self.gas_z_cen - self.gas_x_cen * self.gas_vz_cen
        self.gas_jz_cen = self.gas_vy_cen * self.gas_x_cen - self.gas_y_cen * self.gas_vx_cen
        self.gas_j_mag  = sqrt(self.gas_jx_cen**2. + self.gas_jy_cen**2. + self.gas_jz_cen**2.)






    def write_fits(self):
        print '\tGenerating fits for %s...'%self.aname
        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the momentum measurements in this FITS file."
        prihdr['simname'] = self.simname
        prihdr['scale'] = self.aname.strip('a')
        prihdr['snapfile'] = self.snapfile

        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)

        colhdr = fits.Header()
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_x_cen , self.stars_y_cen , self.stars_z_cen)), header = colhdr , name = 'stars_position'))
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.gas_x_cen , self.gas_y_cen , self.gas_z_cen)), header = colhdr , name = 'gas_position'))
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_jx_cen, self.stars_jy_cen, self.stars_jz_cen)), header = colhdr, name = 'stars_momentum'))
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.gas_jx_cen, self.gas_jy_cen, self.gas_jz_cen)), header = colhdr, name = 'gas_momentum'))


        '''
        master_hdulist.append(fits.ImageHDU(data = self.cube, header = self.cube_hdr, name = 'cam%i_obs_cube'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_int, self.edisp_int]), name = 'cam%i_disp_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.disp_obs, self.edisp_obs]), name = 'cam%i_disp_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_int,self.evel_int]), name = 'cam%i_vel_int'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = array([self.vel_obs,self.evel_obs]), name = 'cam%i_vel_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_obs, name = 'cam%i_ha_obs'%self.camera))
        master_hdulist.append(fits.ImageHDU(data = self.ha_int, name = 'cam%i_ha_int'%self.camera))
        '''

        thdulist = fits.HDUList(master_hdulist)
        print '\tSaving to ' + self.fits_name
        thdulist.writeto(self.fits_name, clobber = True)



        return master_hdulist








def measure_momentum(snapfile, out_sim_dir, nir_cat, nir_disc_cat):
    print 'Measuring momentum for '+ snapfile
    aname = (os.path.basename(snapfile)).split('_')[-1].rstrip('.d')
    simname = snapfile.split('_')[0]
    fits_name = out_sim_dir+'/'+simname+'_'+aname+'_momentum.fits'
    mom = momentum_obj(simname, aname, snapfile, fits_name)
    mom.load()



    in_nir = where(nir_cat[:,0] == aname)[0]
    if len(in_nir) == 0: return
    nir_cat = nir_cat[in_nir[0]]
    nir_disc_cat = nir_disc_cat[in_nir[0]]
    L_disk_s = nir_cat[8:11].astype('float')
    L_disk   = nir_disc_cat[7:10].astype('float')

    mom.load()
    mom.calc_momentum(nir_cat, nir_disc_cat)
    mom.write_fits()








if __name__ == "__main__":

    args = parse()
    import yt

    if args['snap_files'] is not None:
        snaps = [args['snap_files']]
    else:
        snaps = np.asarray(glob.glob("*.d"))

        
    print "Generating Sunrise Input for: ", snaps

    abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(abssnap)

    dirname = os.path.dirname(abssnap)
    simname = os.path.basename(dirname) #assumes directory name for simulation name

    nir_cat_name = simname[0:-2]+'_v2_'+simname[-2:]
    nir_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Nir_simplified_disc_cat.txt', skiprows = 1, dtype='str')
    nir_disc_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Nir_disc_cat.txt', skiprows = 1, dtype='str')


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


    Parallel(n_jobs = -1)(delayed(measure_momentum)(new_snapfiles[i], out_sim_dir, nir_cat, nir_disc_cat) for i in arange(len(new_snapfiles)))



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
















