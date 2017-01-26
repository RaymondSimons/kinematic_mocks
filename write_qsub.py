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


    args = vars(parser.parse_args())
    return args



gal = 'VELA28'



if __name__ == "__main__":
    args = parse()
    if args['gal'] is not None: gal = args['gal']
    else: gal = 'VELA28'

    qsub_direct = '/nobackupp2/rcsimons/momentum_measurements/qsub'



    fname = qsub_direct+'/momentum_'+gal+'.qsub'
    f = open(fname, 'w+')

    f.write('#PBS -S /bin/bash\n')
    f.write('#PBS -l select=1:ncpus=20:model=has\n')
    f.write('#PBS -l walltime=02:00:00\n')
    f.write('#PBS -q normal\n')
    f.write('#PBS -N momentum_measurements\n')
    f.write('#PBS -M rsimons@jhu.edu\n')
    f.write('#PBS -m abe\n')
    f.write('#PBS -o ./out_err/%s.out\n'%gal)
    f.write('#PBS -e ./out_err/%s.err\n'%gal)
    f.write('#PBS -V\n')

    f.write('cd /nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/%s\n'%gal)
    f.write('python /u/rcsimons/scripts/kinematic_mocks/measure_momentum.py > ./out_err/%s.out 2> ./out_err/%s.err\n\n\n'%(gal, gal))

    f.close()
















