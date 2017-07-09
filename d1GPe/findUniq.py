import numpy as np
import pickle
import os
import datetime as dt
storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/forReview/d1GPe'

path1 = storage_home+'/output/'
combs = []
for i in np.arange(0,1693):
	fname = path1+"Combinations_new"+str(0+i)+".pickle"
	fname1 = path1+"Combinations_new_WithConds"+str(0+i)+".pickle"
	if os.path.exists(fname):
		st = os.stat(fname)
		mtime = dt.datetime.fromtimestamp(st.st_mtime)
		if mtime.day >=7 and mtime.month == 7:
			print mtime
			combs.append(pickle.load(open(fname,"r")))
	if os.path.exists(fname1): 
		combs.append(pickle.load(open(fname1,"r")))
seqs = [ np.round(subseq[1],2)  for seq in combs for subseq in seq ]
Allstrs = [ np.array2string(seq) for seq in seqs]
s,uniqueIds = np.unique(Allstrs, return_index=True)
pickle.dump(np.array(seqs)[uniqueIds],open(path1+"uniqueCombsWithConds.pickle","w"))


