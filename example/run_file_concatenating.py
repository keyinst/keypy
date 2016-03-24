import glob, os
from keypy.preprocessing.concatenating_segments import *


##----------- Add info on files you would like to concatenate here ---------------##


##specify folder of files which are to be concatenated
mydir = os.path.join(library_path,"data","input","groups")

##absolute path example
#mydir = os.path.join("E:" + os.path.sep, "data","input","groups")

##specify number of time frames per epoch, number of channels
ntf=512
nch=32

##specify signifier for group participant cond run
#unique identifier
id_start=0
id_end=4

## ending of input files
ending_in = 'asc'

##desired ending of output files
ending_out = 'asci'

##specify signifier for epoch number (currently not used, later version could make sure that epochs are concatenated in particular order)
#ep_start=6
#ep_end=8

#run concatenating
concatenate(mydir, ntf, nch, id_start, id_end, ending_in, ending_out)