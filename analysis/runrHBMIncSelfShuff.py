import numpy as np
import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.fftpack as spfft

import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *
from scipy.stats import ks_2samp


storage_home = os.getcwd() 

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 
pathname = storage_home+'/output/'


labels=["RobustHealthy","AR","TDMild","TDSevere","RobustPathological","SevereAkinesia","RH<->TDMild","TDMild<->TDSevere","RH<->AR","AR<->SA","SA<->RP","InBetween" ]





# Now calculate K-S and generate one script per parameter 

KSPos = np.zeros((6,1)) # K-S statistic
KSNeg = np.zeros((14,1))

MeansClus1Neg = np.zeros((14,1)) # Means
MeansClus2Neg = np.zeros((14,1))
MeansClus1Pos =  np.zeros((6,1))
MeansClus2Pos =  np.zeros((6,1))
VarsClus1Neg = np.zeros((14,1)) # Means
VarsClus2Neg = np.zeros((14,1))
VarsClus1Pos =  np.zeros((6,1))
VarsClus2Pos =  np.zeros((6,1))


terms = ["d1ta","d1ti","d2ta","d2ti","fsita","fsiti","tad2","tata","tati","tid2","tita","titi","stnta","stnti","tistn","tastn","d1ctx","d2ctx","fsictx","stnctx"]
# And their respetive ids in the matrix
ids = [(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,1),(3,3),(3,4),(4,1),(4,3),(4,4),(5,3),(5,4),(4,5),(3,5),(0,7),(1,7),(2,7),(5,7)]
posbin = np.arange(0,19,1)*0.75
negbin = np.arange(-8,0.5,0.3)*0.75


Seeds = pickle.load(open(pathname+"Seeds.pickle","r"))

AllAsPD = pickle.load(open(pathname+"AllAsPDShuff_"+str(Seeds[0])+".pickle","r"))

AllAsH = pickle.load(open(pathname+"AllAsHShuff_"+str(Seeds[0])+".pickle","r"))
clusPD = np.random.randint(0,len(AllAsPD),500)
clusH = np.random.randint(0,len(AllAsH),500)


curr_ind = 0
#for i in xrange(1,22):# Here i represents how many elements to replace by their means
for i in xrange(1,22):# Here i represents how many elements to replace by their means

	for trial,seed1 in enumerate(Seeds):	
		# Calculate 
		params = dict()
		params["clusPD"] = clusPD
		params["clusH"] = clusH
		#params["KSAll"] = KSAll
		params["samplesPD"] = forReplaceInc["samplesPD"].tolist()
		params["samplesH"] = forReplaceInc["samplesH"].tolist()
		#params["tupleSorted"] = forReplaceInc["tupleSorted"]
		params["indicesOrder"] = forReplaceInc["indicesOrder"][:i]
		params["conn"] = forReplaceInc["termsOrder"][:i].tolist()
		params["ind1"] = i
		params["indicesOrdPD"] = forReplaceInc["indicesOrdPD"][trial][:i]
		params["indicesOrdH"] = forReplaceInc["indicesOrdH"][trial][:i]
		params["seed"] = seed1
		sim_name ="replaceHistsByMeansIncSelfShuff_"+str(seed1)
		fh = open(path3 + '%s%d.py'%(sim_name,i),'w')
		fh.write('import sys\n')
		fh.write("sys.path.insert(0,'/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/')\n")
		fh.write('import replaceHistsByMeansIncSelfShuff as rHBM\n')
		fh.write('rHBM.replace(%s)\n'%(params))
		fh.close()


		print sim_name

		fh = open(path3 + '%s.jdf'%(sim_name),'w')
		content = [
		'#!/bin/bash\n',
		#   '#PBS -t 0-%d\n'%(len(comb_pars)-1),
		#'#SBATCH -o /home/j.bahuguna/homology/vAModel/output/test_job_.out',
		'#SBATCH --output=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/output/rHBM_%s%d.out'%(sim_name,i),
		'#SBATCH --error=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/output/rHBM_%s%d.err'%(sim_name,i),
		#'#PBS -j oe\n',
		'#SBATCH --job-name=%s\n'%(sim_name),
		'#SBATCH --mem-per-cpu=3500mb\n',
		'#SBATCH --time 5:00:00\n',
		#'#SBATCH -p long \n',
		#'#SBATCH --output=%s%s_%s.txt\n'%(path3,sim_name,postfix),
		#'export PYTHONPATH=/clustersw/cns/nest-2.2.2/lib/python2.7/dist-packages/:$PYTHONPATH\n',
		'python /home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/jdf/%s%d.py\n'%(sim_name,i),
		#'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
		]
		fh.writelines(content)
		fh.close()
		filename = path3 + '%s.jdf'%(sim_name)	
		os.system('sbatch  %s'%filename )

