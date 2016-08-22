import os
import glob


mcrx = glob.glob("mcrx.fits")



print os.path.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))

basepath = os.path.dirname(os.getcwd().replace('/ifu','').replace('_sunrise',''))
print basepath[len(basepath)-5::]






#print os.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
