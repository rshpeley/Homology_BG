import numpy as np
import itertools
import shutil
import os
import pickle

sim_name = "PDFH_slow"
#storage_home = os.environ.get('STORAGE_HOME')
storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy'

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 

input_range = np.arange(0,9,1)
Allcombs = []
pathname = storage_home+'/output/'
allnames = os.listdir(pathname)
indices = np.random.random_integers(0,len(allnames)-1,1100)
#indices = pickle.load(open("RandomIndicesGood.pickle","r"))
#for filename in np.array(allnames)[indices][::1]:
#for no in np.arange(0,3000):
for no in np.arange(0,7500):
	#filename = "Allcombs_Comb_NewSig"+str(no)+"_.pickle"
	filename = "Allcombs_Tight_withThaComb_"+str(no)+"_.pickle"
#for filename in np.array(allnames):
	if os.path.isfile(pathname+filename) == True:
		temp = pickle.load(open(pathname+filename,"r"))
		if len(temp) > 0:
			for x in temp:
				#if len(x) == 17:
				if len(x) == 20:
					Allcombs.append(np.array(x))

pickle.dump(Allcombs,open(path5+"AllcombsQuick.pickle","w"))
sim_name = "GenABs_quick_"
fh = open(path3 + '%s.py'%(sim_name),'w')
fh.write('import sys\n')
fh.write("sys.path.insert(0,'/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/')\n")
fh.write('import genValidABs as gVAB\n')
fh.write('gVAB.run()\n')
fh.close()

print sim_name

fh = open(path3 + '%s.jdf'%(sim_name),'w')
content = [
'#!/bin/bash\n',
#   '#PBS -t 0-%d\n'%(len(comb_pars)-1),
#'#SBATCH -o /home/j.bahuguna/homology/vAModel/output/test_job_.out',
'#SBATCH --output=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/output/gVAB_%s.out'%(sim_name),
'#SBATCH --error==/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/output/gVAB_%s.err'%(sim_name),
#'#PBS -j oe\n',
'#SBATCH --job-name=%s\n'%(sim_name),
'#SBATCH --mem-per-cpu=3500mb\n',
'#SBATCH --time 24:00:00\n',
#'#SBATCH -p long \n',
#'#SBATCH --output=%s%s_%s.txt\n'%(path3,sim_name,postfix),
#'export PYTHONPATH=/clustersw/cns/nest-2.2.2/lib/python2.7/dist-packages/:$PYTHONPATH\n',
'python /home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/jdf/%s.py\n'%(sim_name),
#'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
]
fh.writelines(content)
fh.close()
filename = path3 + '%s.jdf'%(sim_name)	
os.system('sbatch  %s'%filename )

