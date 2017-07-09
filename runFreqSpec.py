import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
#mpl.use('Agg')
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.fftpack as spfft

import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *
#import tree_anal
#import paramsearchGA_DopDep as psGA
import paramsearchGA_DopDep_nonlinear as psGA
import knownUnknownParams as p
#from sklearn.cluster import KMeans
import funcs as fs 

knownparams = p.params["known"]


terms = ["d1ta","d1ti","d2ta","d2ti","fsita","fsiti","tad2","tata","tati","tid2","tita","titi","stnta","stnti","tistn","tastn","jc1","jc2","jfsictx","jstnctx"]
terms1 = ["D1","D2","FSI","TA","TI","STN","GPi","Ctx"]



storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams'

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 


clusterIndsPD = pickle.load(open(path5+"clustersPDMore5.pickle","r"))
for i,clusPD in enumerate(clusterIndsPD):
	if len(clusPD) > 0 :
		print i

		params=dict()
		params["clusPD"] = clusPD	
		params["i"] = i
		sim_name = "freqSpec_"+str(i)
		
		fh = open(path3 + '%s.py'%(sim_name),'w')
		fh.write('import sys\n')
		fh.write("sys.path.insert(0,'/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/')\n")
		fh.write('import analysisFreqSpec as aFS\n')
		fh.write('aFS.analFreqSpec(%s)\n'%(params))
		fh.close()


		fh = open(path3 + '%s.jdf'%(sim_name),'w')
		content = [
		'#!/bin/bash\n',
		#   '#PBS -t 0-%d\n'%(len(comb_pars)-1),
		#'#SBATCH -o /home/j.bahuguna/homology/vAModel/output/test_job_.out',
		'#SBATCH --output=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/output/FS_%s.out'%(sim_name),
		'#SBATCH --error=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/output/FS_%s.err'%(sim_name),
		#'#PBS -j oe\n',
		'#SBATCH --job-name=%s\n'%(sim_name),
		'#SBATCH --mem-per-cpu=3500mb\n',
		'#SBATCH --time 24:00:00\n',
		#'#SBATCH -p long \n',
		#'#SBATCH --output=%s%s_%s.txt\n'%(path3,sim_name,postfix),
		#'export PYTHONPATH=/clustersw/cns/nest-2.2.2/lib/python2.7/dist-packages/:$PYTHONPATH\n',
		'python /home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/jdf/%s.py\n'%(sim_name),
		#'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
		]
		fh.writelines(content)
		fh.close()
		filename = path3 + '%s.jdf'%(sim_name)	
		os.system('sbatch  %s'%filename )

