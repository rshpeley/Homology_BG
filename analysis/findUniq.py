import numpy as np
import pickle
import os

storage_home = os.getcwd()"../PD/output/" # Or healthy 

path1 = storage_home+'/output/'
combs = []
for i in np.arange(0,409):
	fname = path1+"Combinations_new"+str(0+i)+".pickle"
	if os.path.exists(fname):
		combs.append(pickle.load(open(fname,"r")))

seqs = [ np.round(subseq[1],2)  for seq in combs for subseq in seq ]
Allstrs = [ np.array2string(seq) for seq in seqs]
s,uniqueIds = np.unique(Allstrs, return_index=True)
pickle.dump(np.array(seqs)[uniqueIds],open(path1+"uniqueCombs.pickle","w"))


