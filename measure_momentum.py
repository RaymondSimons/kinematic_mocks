import numpy as np
from numpy import *
import os, sys, argparse
import glob
import astropy
from astropy.io import fits

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
        
        print "Sunrise directory: ", snap_dir
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


    # Loop over snapshot to generate cameras and projection plots, 
    # parallelization happens while generating the plots.
    galprops_file = simname+'_galprops.npy'
    galprops = np.load(galprops_file)[()]



    for snapfile in new_snapfiles:

        aname = (os.path.basename(snapfile)).split('_')[-1].rstrip('.d')

        print "Timestep name: ", aname

        snap_dir = os.path.dirname(snapfile.rstrip('_sunrise')) #os.path.join(simname+'_'+aname+'_sunrise')

        print "Sunrise directory: ", snap_dir
        assert os.path.lexists(snap_dir)



        out_sim_dir = os.path.join('/nobackupp2/rcsimons/momentum_measurements/', simname)
        print os.path.lexists(out_sim_dir)
        if not os.path.lexists(out_sim_dir):
            os.mkdir(out_sim_dir)                    

        out_dir = os.path.join(out_sim_dir, snap_dir)
        print os.path.lexists(out_dir)
        if not os.path.lexists(out_dir):
            os.mkdir(out_dir)


        if os.path.abspath(snapfile) not in galprops['snap_files']: continue
        idx = np.argwhere(galprops['snap_files']==os.path.abspath(snapfile))[0][0]





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
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Gas circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Gas circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Gas circularity_r'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Young Stars circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Young Stars circularity_r'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Old Stars circularity_z'))
        col_list.append(fits.ImageHDU(data = np.empty((10,10)), name = 'Old Stars circularity_r'))

        thdulist = fits.HDUList(col_list)
        thdulist.writeto(out_dir+'/'+snap_dir+'_kinematics.fits', clobber = True)




