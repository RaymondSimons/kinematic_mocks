import os
import glob


mcrx = glob.glob("mcrx.fits")



print os.path.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
print os.path.dirname(os.getcwd().replace('/ifu','').replace('_sunrise',''))





#print os.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
