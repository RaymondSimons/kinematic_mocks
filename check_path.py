import os
import glob


mcrx = glob.glob("mcrx.fits")



kmap_dir = '/nobackupp2/rcsimons/sunrise_testing/kmaps'


basepath = os.path.dirname(os.getcwd().replace('/ifu','').replace('_sunrise',''))



simname = basepath[len(basepath)-6::]
snapname = os.path.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))

print simname, snapname






#print os.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
