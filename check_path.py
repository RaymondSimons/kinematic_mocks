import os
import glob


mcrx = glob.glob("mcrx.fits")

print os.getcwd().replace('/ifu','').replace('_sunrise','')
#print os.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
