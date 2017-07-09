import numpy as np
import pickle
import os

storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/forReview/stnta_zero'

path1 = storage_home+'/output/'
combs = []
for i in np.arange(290,409):
	fname = path1+"Combinations_new"+str(0+i)+".pickle"
	if os.path.exists(fname):
		combs.append(pickle.load(open(fname,"r")))

seqs = [ np.round(subseq[1],2)  for seq in combs for subseq in seq ]
Allstrs = [ np.array2string(seq) for seq in seqs]
s,uniqueIds = np.unique(Allstrs, return_index=True)
pickle.dump(np.array(seqs)[uniqueIds],open(path1+"uniqueCombs.pickle","w"))


