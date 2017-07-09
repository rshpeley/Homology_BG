import sys
import os
import numpy as np
start = sys.argv[1]
stop = sys.argv[2]
print start
print stop
for id in np.arange(int(start),int(stop),1):
	cmd = "scancel "+str(id)
	os.system(cmd)


