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
    fsh = open(qsub_direct+'/'+gal+'_submit_momentum.sh', 'w+')

    snaps = np.asarray(glob.glob("/nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/%s/*.d"%gal))

    for sn in snaps:
        dname = os.path.basename(sn)
        aname = dname.split('_')[2].strip('.d')
        fname = qsub_direct+'/momentum_%s_%s.qsub'%(gal, aname)
        fsh.write('sh '+fname+'\n')
        f = open(fname, 'w+')
        f.write('#PBS -S /bin/bash\n')
        f.write('#PBS -l select=1:ncpus=24:model=has\n')
        f.write('#PBS -l walltime=01:00:00\n')
        f.write('#PBS -q normal\n')
        f.write('#PBS -N %s_%s_momentum\n'%(gal, aname))
        f.write('#PBS -M rsimons@jhu.edu\n')
        f.write('#PBS -m abe\n')
        f.write('#PBS -o ./out_err/%s_%s_pbs.out\n'%(gal, aname))
        f.write('#PBS -e ./out_err/%s_%s_pbs.err\n'%(gal, aname))
        f.write('#PBS -V\n')

        f.write('cd /nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/%s\n'%gal)
        f.write('python /u/rcsimons/scripts/kinematic_mocks/measure_momentum.py %s> /nobackupp2/rcsimons/momentum_measurements/qsub/out_err/%s_%s.out 2> /nobackupp2/rcsimons/momentum_measurements/qsub/out_err/%s_%s.err\n\n\n'%(sn, gal, aname, gal, aname))

        f.close()

    fsh.close()






















