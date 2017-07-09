import numpy as np
import itertools
import shutil
import os
import pickle
import sys
import time
sim_name = "Dist_PD_slow"
storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy'
path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/jdf/' # Path to jdf files
path5 = storage_home+'/output/' # Path for output 
combs = []

for i in np.arange(0,14900,1):
	fname = path5+"Allcombs_Tight_withThaComb_"+str(0+i)+"_.pickle"
	if os.path.exists(fname):
		temp = pickle.load(open(fname,"r"))
		if len(temp) > 0:
			for x in temp:
				combs.append(x)


pickle.dump(combs,open(path5+"uniqueCombs.pickle","w"))

