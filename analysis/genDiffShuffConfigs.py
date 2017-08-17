import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
#mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.fftpack as spfft
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import scipy.stats as sp1
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *

from mpl_toolkits.axes_grid.inset_locator import inset_axes,InsetPosition
# Read all the 3 clusters
storage_home = os.getcwd() 

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 
pathname = storage_home+'/output/'

AllAsPD = pickle.load(open(storage_home+"../PD/output/AllAs.pickle","r"))
AllBsPD = pickle.load(open(storage_home+"../PD/output/AllBs.pickle","r"))

AllAsH = pickle.load(open(storage_home+"../healthy/output/AllAs.pickle","r"))
AllBsH = pickle.load(open(storage_home+"../healthy/output/AllBs.pickle","r"))

# names of the free parameters .. 
terms = ["d1ta","d1ti","d2ta","d2ti","fsita","fsiti","tad2","tata","tati","tid2","tita","titi","stnta","stnti","tistn","tastn","d1ctx","d2ctx","fsictx","stnctx"]
# And their respetive ids in the matrix
ids = [(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,1),(3,3),(3,4),(4,1),(4,3),(4,4),(5,3),(5,4),(4,5),(3,5),(0,7),(1,7),(2,7),(5,7)]
posbin = np.arange(0,19,0.5)
negbin = np.arange(-8,0.6,0.5)


AllAsPDShuff=[]
AllBsPDShuff=[]
AllAsHShuff=[]
AllBsHShuff=[]

Seeds = np.random.randint(0,99999990,10)
pickle.dump(Seeds,open(pathname+"Seeds.pickle","w"))

# For different trials - different seed
for x in xrange(10): # 10 trials for now
	seed = Seeds[x]
	print "seed",seed
	np.random.seed(seed)	
	AsPDShuff=np.copy((AllAsPD))
	BsPDShuff=np.copy((AllBsPD))
	AsHShuff=np.copy((AllAsH))
	BsHShuff=np.copy((AllBsH))

	for i in xrange(len(terms)):
		if i < 16:
			ser12 = np.array(AllAsPD)[:,ids[i][0],ids[i][1]]
			print "unshuff mean",np.mean(np.array(AllAsPD)[:,ids[i][0],ids[i][1]] )
			np.random.shuffle(ser12)
			
			AsPDShuff[:,ids[i][0],ids[i][1]] = ser12
			ser1 = np.array(AllAsPD)[:,ids[i][0],ids[i][1]] # Can be commented once you are sure, things are working correctly
			print "shuff mean",np.mean(AsPDShuff[:,ids[i][0],ids[i][1]]) # Can be commented

			ser12 = np.array(AllAsH)[:,ids[i][0],ids[i][1]]
			print "H -unshuff mean",np.mean(ser12)
			np.random.shuffle(ser12)
			
			AsHShuff[:,ids[i][0],ids[i][1]] = ser12
			print "H -shuff mean",np.mean(AsHShuff[:,ids[i][0],ids[i][1]]) # Can be commented



		else:
			ser12 = np.array(AllBsPD)[:,0,ids[i][0]]
			print "unshuff mean",np.mean(ser12)
			np.random.shuffle(ser12)
			BsPDShuff[:,0,ids[i][0]] = ser12
			ser1 = np.array(AllBsPD)[:,0,ids[i][0]]
			print "shuff mean",np.mean(BsPDShuff[:,0,ids[i][0]]) # Can be commented

			ser12 = np.array(AllBsH)[:,0,ids[i][0]]
			print "H -unshuff mean",np.mean(ser12)

			np.random.shuffle(ser12)
			BsHShuff[:,0,ids[i][0]] = ser12
			print "H - shuff mean",np.mean(BsHShuff[:,0,ids[i][0]]) # Can be commented


	
	pickle.dump(AsPDShuff,open(pathname+"AllAsPDShuff_"+str(seed)+".pickle","w"))
	pickle.dump(AsHShuff,open(pathname1+"AllAsHShuff_"+str(seed)+".pickle","w"))
	pickle.dump(BsPDShuff,open(pathname+"AllBsPDShuff_"+str(seed)+".pickle","w"))
	pickle.dump(BsHShuff,open(pathname1+"AllBsHShuff_"+str(seed)+".pickle","w"))
	

