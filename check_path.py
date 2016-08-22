import os
import glob


mcrx = glob.glob("mcrx.fits")



kmap_dir = '/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2'


basepath = os.path.dirname(os.getcwd().replace('/ifu','').replace('_sunrise',''))



simname = basepath[len(basepath)-6::]
snapname = os.path.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))

print simname, snapname

kmap_name = '/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname+'/'+snapname
print kmap_name
print os.path.lexists('/nobackupp2/rcsimons/sunrise_testing/kmaps/VELA_v2/'+simname)






#print os.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
