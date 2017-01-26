import numpy as np
import glob
import os, sys, argparse



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

    parser.add_argument('gal', nargs='?', default=None, help='Snapshot files to be analyzed.')

    #parser.add_argument('-s', '--snap_base', default='10MpcBox_csf512_',
    #                    help='Base of the snapshots file names.') 

    #parser.add_argument('-d', '--distance', default=100000, type=float,
    #                    help='Distance between cameras and the center of the galaxy (in [kpc]).')

    #parser.add_argument('--no_export',action='store_true',
    #                    help='Do not export data to fits for Sunrise.') 

    args = vars(parser.parse_args())
    return args



gal = 'VELA28'



if __name__ == "__main__":
    args = parse()

    if args['gal'] is not None: gal = args['gal']
    else: gal = 'VELA28'

    snaps = glob.glob('/nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/'+gal+'/*.d')
    print snaps
    for sn in snaps:
        aname = sn.split('_')[-1].rstrip('.d').strip('a')
        print aname






#PBS -S /bin/bash
#PBS -l select=1:ncpus=1:model=has
#PBS -l walltime=04:00:00
#PBS -q normal
#PBS -N sunrise_export
#PBS -M gsnyder@stsci.edu
#PBS -m abe
#PBS -o sunrise_export_a0.200pbs.out
#PBS -e sunrise_export_a0.200pbs.err
#PBS -V
