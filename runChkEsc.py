'''
import pickle
import paramsearchGA_DopDep as psGA
path = "/storage/users/jb1022/output/"
delay = 1.0
elmt = psGA.paramSearch(delay)
pickle.dump(elmt,open(path+"Combinations_DopDep.pickle","w"))
# Running on the non-grid machine
'''
import numpy as np
import itertools
import shutil
import os
import pickle

sim_name = "check_Esc_PD"
#storage_home = os.environ.get('STORAGE_HOME')
storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams'

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 

ctxinp = np.arange(0,15,1)

for ip in ctxinp:
	params = dict()
	params["ipAmp"] = ip
	sim_name ="check_Esc_PD"+str(ip)
	fh = open(path3 + '%s.py'%(sim_name),'w')
	fh.write('import sys\n')
	fh.write("sys.path.insert(0,'/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/')\n")
	fh.write('import checkEscande as cEs\n')
	fh.write('cEs.calcFI(%s)\n'%(params))
	fh.close()
	
	print sim_name

	fh = open(path3 + '%s.jdf'%(sim_name),'w')
	content = [
	'#!/bin/bash\n',
	#   '#PBS -t 0-%d\n'%(len(comb_pars)-1),
	#'#SBATCH -o /home/j.bahuguna/homology/vAModel/output/test_job_.out',
	'#SBATCH --output=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/output/chkEs_%s.out'%(sim_name),
	'#SBATCH --error=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/output/chkEs_%s.err'%(sim_name),
	#'#PBS -j oe\n',
	'#SBATCH --job-name=%s\n'%(sim_name),
	'#SBATCH --mem-per-cpu=3500mb\n',
	'#SBATCH --time 8:00:00\n',
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

