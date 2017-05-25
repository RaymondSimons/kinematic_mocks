import pyfits
import tarfile
import glob
from glob import glob


for i in arange(1,35):
    if i!=18:
        f = open('/nobackupp2/rcsimons/catalogs/image_cats/VELA%.2i_image.cat'%(i))
        path = '/nobackupp2/gfsnyder/VELA_sunrise/Runs/VELA_v2/VELA%.2i'%(i)
        fls = glob(path+'/*sunrise.tar')
        for fl in fls:
            print fl
