# This is the main script that generates scripts with all combinations of the parameters 

import numpy as np
import itertools
import shutil
import os
import pickle
import time

sim_name = "Dist_H_slow"
storage_home = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy"
path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/jdf/' # Path to jdf files
path5 = storage_home+'/output/' # Path for output 
combs = []
#combs.append(pickle.load(open("Combfiles/Combinations_10.pickle","r")))
for i in np.arange(0,516,2):
	fname = path5+"Combinations_new"+str(0+i)+".pickle"
	if os.path.exists(fname):
		combs.append(pickle.load(open(fname,"r")))

#combs = np.array(combs)
seqs = [ subseq[1]  for seq in combs for subseq in seq ]
Allstrs = [ np.array2string(seq) for seq in seqs]
s,uniqueIds = np.unique(Allstrs, return_index=True)
pickle.dump([seqs,uniqueIds],open("uniqueCombs1.pickle","w"))

Ucombs = pickle.load(open("uniqueCombs1.pickle","r"))
seqs = Ucombs[0]
uniqueIds = Ucombs[1]

batch_pars = {
    'ind' :np.array(seqs)[uniqueIds] # These are all the unique combinations of indices in range of 16 parameters ,for which criteria was met
}
pickle.dump(batch_pars['ind'],open(path5+"uniqueLastIter.pickle","w"))

Samples = np.arange(0,len(batch_pars['ind']),1)
for x in np.arange(0,len(batch_pars['ind']),150):

	#offset =4600
	offset =x

	randSamps = Samples[offset:offset+150] # And then next round, next 1000
	for jj,ii in enumerate(randSamps):
		postfix = "Comb_"+str(jj+offset)+"_"
		#postfix = "Comb_NewSig"+str(jj+0)+"_"
		fh = open(path1 + '%s_%s.py'%(sim_name,postfix),'w')
		fh.write('import sys\n')
		fh.write("sys.path.insert(0,'/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/')\n")
		fh.write('import findDist as fD\n')
		fh.write('prefix = "%s"\n'%postfix)
		#fh.write('prefix = "sim_1_%d_%d"\n'%(g,eta))
		fh.write('pars = dict()\n')
		fh.write('pars["%s"] = ['%('ind')) 
		fh.write(",".join(["%d" %i for i in batch_pars['ind'][ii]]))
		fh.write("]\n")	
		fh.write('fD.run(prefix,pars)\n')
		#fh.write('net1.sim()\n')
		fh.close()
		print sim_name,'%s'%(postfix)

	'''
	fname = path1 + '%s.log'%sim_name
	file(fname,'w').write(str(batch_pars))
	# copy generic params_d file and store with modified name;
	# this will be later used by the batch scripts
	shutil.copy(path3,path4)
	'''


	#create jdf file
	for jj,ii in enumerate(randSamps):
	#fh = open(path2 + '%s.jdf'%sim_name,'w')
		postfix = "Comb_"+str(jj+offset)+"_"
		#postfix = "Comb_NewSig"+str(jj+0)+"_"
		fh = open(path2 + '%s.jdf'%(sim_name),'w')
		content = [
		'#!/bin/bash\n',
		#   '#PBS -t 0-%d\n'%(len(comb_pars)-1),
		#'#SBATCH -o /home/j.bahuguna/homology/vAModel/output/test_job_.out',
		'#SBATCH --output=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/output/PDF_%s.out'%(sim_name),
		'#SBATCH --error=/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/output/PDF_%s.err'%(sim_name),
		#'#PBS -j oe\n',
		'#SBATCH --job-name=%s\n'%(sim_name),
		'#SBATCH --mem-per-cpu=3500mb\n',
		'#SBATCH --time 1:00:00\n',
		#'#SBATCH -p long \n',
		#'#SBATCH --output=%s%s_%s.txt\n'%(path3,sim_name,postfix),
		#'export PYTHONPATH=/clustersw/cns/nest-2.2.2/lib/python2.7/dist-packages/:$PYTHONPATH\n',
		'python /home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/scripts/%s_%s.py\n'%(sim_name,postfix),
		#'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
		]
		fh.writelines(content)
		fh.close()
		filename = path2 + '%s.jdf'%(sim_name)	
		os.system('sbatch  %s'%filename )

	# After posting 100 jobs, sleep for an hour and then post again 100 jobs
	time.sleep(60*7)		

	

