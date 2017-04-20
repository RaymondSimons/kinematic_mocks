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
    if args['gal'] is not None: 
        gal = args['gal']

        qsub_direct = '/nobackupp2/rcsimons/data/kin_maps/%s/qsub'%gal
        fsh = open('/nobackupp2/rcsimons/runs_files/'+gal+'_submit_mockcubes.sh', 'w+')

        fname = qsub_direct+'/%s_mockcubes.qsub'%gal
        fsh.write('qsub '+fname+'\n')
        f = open(fname, 'w+')
        f.write('#PBS -S /bin/bash\n')
        f.write('#PBS -l select=1:ncpus=24:model=has\n')
        f.write('#PBS -l walltime=01:30:00\n')
        f.write('#PBS -q normal\n')
        f.write('#PBS -N %s_kmap\n'%gal)
        f.write('#PBS -M rsimons@jhu.edu\n')
        f.write('#PBS -m abe\n')
        f.write('#PBS -o %s/%s_kmap_pbs.out\n'%(qsub_direct, gal))
        f.write('#PBS -e %s/%s_kmap_pbs.err\n'%(qsub_direct, gal))
        f.write('#PBS -V\n')

        f.write('cd /nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/%s\n'%gal)

        comm_1 = 'python /u/rcsimons/scripts/kinematic_mocks/make_kin_fits.py'
        outf   = '%s/%s_kmap.out'%(qsub_direct, gal)
        errf   = '%s/%s_kmap.err'%(qsub_direct, gal)

        f.write('%s > %s 2> %s \n'%(comm_1, outf, errf))

        f.write('tar -zcvf /nobackupp2/rcsimons/data/kin_maps/%s_kmaps.tar.gz /nobackupp2/rcsimons/data/kin_maps/%s/*_kmap.fits\n\n'%(gal, gal))

        f.close()

        fsh.close()
    else:
        print 'No galaxy name given'
