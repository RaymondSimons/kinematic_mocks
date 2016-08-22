import os
import glob


mcrx = glob.glob("mcrx.fits")

print os.basename(os.getcwd().replace('/ifu','').replace('_sunrise',''))
