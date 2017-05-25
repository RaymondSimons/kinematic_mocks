import pyfits
import tarfile
import glob
from glob import glob


for i in arange(1,2):
    if i!=18:
        f = open('/nobackupp2/rcsimons/catalogs/image_cats/VELA%.2i_image.cat'%(i),'w+')
        path = '/nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/VELA%.2i'%(i)
        print path
        fls = glob(path+'/*/images/*sunrise.tar')
        fls.sort()
        for fl in fls:
            s_strt = fl.find('_sunrise.tar')-3
            s_stop = fl.find('_sunrise.tar')
            scale = fl[s_strt:s_stop]
            f.write(scale + '\t')
            with tarfile.open(fl) as tf:
                for SB in ['SB25', 'SB27', 'SB00']:
                    for inst in ['ACS-F606W','ACS-F775W', 'ACS-F850W', 'WFC3-F105W', 'WFC3-F125W', 'WFC3-F160W']:
                        for cam_n in arange(19):
                            fits_name = 'images_VELA%.2i_a0.%s_sunrise/VELA%.2i_a0.%s_sunrise_cam%.2i_%s_%s.fits'%(i, scale, i, scale, cam_n, inst, SB)
                            print fits_name
                            f.write(SB+'\t'+inst+'\t'+'%.2i'%cam_n+'\t')
                            #fits_file = tf.getmember(tf)    for entry in tf:
            f.write('\n')



        f.close()
