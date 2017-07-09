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

sim_name = "GA_Healthy_Mallet_wojc1jc2"
#storage_home = os.environ.get('STORAGE_HOME')
storage_home = '/users/bahuguna/2ndManuscript/20paramsfree/checkRelSWAACT/jc1NoEqjc2/test/withOsc/Mallet/strictHealthy/wojc1jc2'

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 

postfix =""
fh = open(path3 + '%s_%s.py'%(sim_name,postfix),'w')
fh.write('import sys\n')
fh.write("sys.path.insert(0,'/users/bahuguna/2ndManuscript/20paramsfree/checkRelSWAACT/jc1NoEqjc2/test/withOsc/Mallet/strictHealthy/wojc1jc2/')\n")
fh.write('import paramsearchGA_DopDep_nonlinear as psGA\n')
fh.write('delay = 1.0\n')
fh.write('psGA.paramSearch(delay)\n')
fh.close()
print sim_name

fh = open(path3 + '%s_%s.jdf'%(sim_name,postfix),'w')
content = [
'#!/bin/bash\n',
#'#SBATCH -o %s\n',
'#PBS -j oe\n',
    '#PBS -d %s\n'%path3,
    '#PBS -j oe\n',
    '#PBS -N %s_%s\n'%(sim_name,postfix),
    '#PBS -l mem=3500mb\n',
    '#PBS -l walltime=24:00:00\n',
    '#PBS -l procs=4\n'
##'export PYTHONPATH=/clustersw/cns/nest-2.2.2/lib/python2.7/dist-packages/:$PYTHONPATH\n',
'python /users/bahuguna/2ndManuscript/20paramsfree/checkRelSWAACT/jc1NoEqjc2/test/withOsc/Mallet/strictHealthy/wojc1jc2/jdf/%s_%s.py\n'%(sim_name,postfix),
#'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
]
fh.writelines(content)
fh.close()
filename = path3 + '%s_%s.jdf'%(sim_name,postfix)	
os.system('qsub  %s'%filename )

