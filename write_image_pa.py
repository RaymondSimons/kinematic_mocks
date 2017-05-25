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
            f.write(scale+'\n')
        f.close()
