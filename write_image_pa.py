import pyfits
import tarfile
import glob
from glob import glob


for i in arange(17, 26):
    print 'VELA%.2i'%(i)
    if i!=18:
        f = open('/nobackupp2/rcsimons/catalogs/image_cats/VELA%.2i_image.cat'%(i),'w+')
        f.write('#(0) scale\n#(1) surface brightness\n#(2) instrument\n#(3) PA camera 0\n#(4) semi-minor axis camera 0\n#(5) semi-major axis camera 0\n#(6) PA camera 1\n#(7) semi-minor axis camera 1\n#(8) semi-major axis camera 1\n#(9:) etc. \n\n\n\n')
        path = '/nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/VELA%.2i'%(i)
        print path
        fls = glob(path+'/*/images/*sunrise.tar')
        fls.sort()
        for fl in fls:
            print '\t', fl
            s_strt = fl.find('_sunrise.tar')-3
            s_stop = fl.find('_sunrise.tar')
            scale = fl[s_strt:s_stop]
            with tarfile.open(fl) as tf:
                for SB in ['SB25', 'SB27']:
                    for inst in ['ACS-F606W','ACS-F775W', 'ACS-F850LP', 'WFC3-F105W', 'WFC3-F125W', 'WFC3-F160W']:
                        f.write(scale + '\t'+ SB+'\t'+inst+'\t')
                        cams = 19
                        pa = zeros(cams)*nan
                        semiminor = zeros(cams)*nan
                        semimajor = zeros(cams)*nan
                        for cam_n in arange(cams):
                            try:
                                fits_name = 'images_VELA%.2i_a0.%s_sunrise/VELA%.2i_a0.%s_sunrise_cam%.2i_%s_%s.fits'%(i, scale, i, scale, cam_n, inst, SB)
                                fits_file = tf.getmember(fits_name)
                                file_obj = tf.extractfile(fits_file)
                                data = pyfits.open(file_obj)
                                pa[cam_n] = data['LotzMorphMeasurements'].header['ORIENT']*180/pi
                                semiminor[cam_n] = data['PhotUtilsMeasurements'].header['SMINSIG']
                                semimajor[cam_n] = data['PhotUtilsMeasurements'].header['SMAJSIG']
                            except:
                                pass
                                #print 'No measurements found for %s %s %s'%(SB, inst, cam_n)
                            f.write('%.2f'%pa[cam_n]+'\t'+'%.3f'%semiminor[cam_n]+'\t'+'%.3f'%semimajor[cam_n]+'\t')
                        f.write('\n')
                                #fits_file = tf.getmember(tf)    for entry in tf:



        f.close()
