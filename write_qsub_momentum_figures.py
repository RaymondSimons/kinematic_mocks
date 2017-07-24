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


if __name__ == "__main__":
    args = parse()
    if args['gal'] is not None: gal = args['gal']
    else: gal = 'VELA28'

    qsub_direct = '/nobackupp2/rcsimons/momentum_measurements/qsub'


    fsh = open(qsub_direct+'/submit_momentum_figures.sh', 'w+')
    gals = ['VELA01', 'VELA06', 'VELA07', 'VELA11', 'VELA15', 'VELA17', 'VELA20',
            'VELA21', 'VELA24', 'VELA27', 'VELA28', 'VELA29', 'VELA33', 'VELA34']
    gals = ['VELA%.2i'%(i+1) for i in np.arange(35)]

    for gal in gals:
        print gal
        fname = qsub_direct+'/momentum_figures_%s.qsub'%gal
        fsh.write('qsub '+fname+'\n')
        f = open(fname, 'w+')
        f.write('#PBS -S /bin/bash\n')
        f.write('#PBS -l select=1:ncpus=24:model=has\n')
        f.write('#PBS -l walltime=00:45:00\n')
        f.write('#PBS -q normal\n')
        f.write('#PBS -N %s_mom_figures\n'%gal)
        f.write('#PBS -M rsimons@jhu.edu\n')
        f.write('#PBS -m abe\n')
        f.write('#PBS -o ./out_err/%s_momentum_figures_pbs.out\n'%gal)
        f.write('#PBS -e ./out_err/%s_momentum_figures_pbs.err\n'%gal)
        f.write('#PBS -V\n')

        f.write('python /u/rcsimons/scripts/kinematic_mocks/momentum_heatmaps.py %s> /nobackupp2/rcsimons/momentum_measurements/qsub/out_err/%s_momentum_figures.out 2> /nobackupp2/rcsimons/momentum_measurements/qsub/out_err/%s_momentum_figures.err\n\n\n'%(gal, gal, gal))

        f.close()

    fsh.close()
