
import numpy as np
import itertools
import shutil
import os
import pickle

sim_name = "PDFPD_test"
#storage_home = os.environ.get('STORAGE_HOME')
storage_home = os.getcwd()+"../PD" # Or healthy

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 

input_range = np.arange(0,9,1.)
for ip in input_range[:1]:
	params = dict()
	params["ipAmp"] = ip
	sim_name = "PDFPD_with20Hz_"+str(ip)
	fh = open(path3 + '%s.py'%(sim_name),'w')
	fh.write('import sys\n')
	pathcmd ='sys.path.insert(0,'+"'"+storage_home+"'"+')\n'
	fh.write(pathcmd)

	fh.write('import checkPDFeaturesStrRed as cPDF\n')
	fh.write('cPDF.classification(%s)\n'%(params))
	fh.close()
	
	print sim_name

	fh = open(path3 + '%s.jdf'%(sim_name),'w')
	content = [
	'#!/bin/bash\n',
	#   '#PBS -t 0-%d\n'%(len(comb_pars)-1),
	#'#SBATCH -o /home/j.bahuguna/homology/vAModel/output/test_job_.out',
	'#SBATCH --output=%s/output/PDF_%s.out'%(storage_home,sim_name),
	'#SBATCH --error=%s/output/PDF_%s.err'%(storage_home,sim_name),
	#'#PBS -j oe\n',
	'#SBATCH --job-name=%s\n'%(sim_name),
	'#SBATCH --mem-per-cpu=3500mb\n',
	'#SBATCH --time 24:00:00\n',
	#'#SBATCH -p long \n',
	#'#SBATCH --output=%s%s_%s.txt\n'%(path3,sim_name,postfix),
	#'export PYTHONPATH=/clustersw/cns/nest-2.2.2/lib/python2.7/dist-packages/:$PYTHONPATH\n',
	'python %s/jdf/%s.py\n'%(storage_home,sim_name),
	#'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
	]
	fh.writelines(content)
	fh.close()
	filename = path3 + '%s.jdf'%(sim_name)	
	os.system('sbatch  %s'%filename )

