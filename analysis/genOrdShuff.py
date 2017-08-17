import pickle
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
#mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.fftpack as spfft

import matplotlib.cm as cm
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *
#import paramsearchGA_DopDep as psGA
import paramsearchGA_DopDep_nonlinear_moreTight as psGA
from sklearn.cluster import KMeans
import knownUnknownParamsStrRed as p

terms = ["d1ta","d1ti","d2ta","d2ti","fsita","fsiti","tad2","tata","tati","tid2","tita","titi","stnta","stnti","tistn","tastn","d1ctx","d2ctx","fsictx","stnctx"]
samplesPD = np.zeros((len(terms),2 )) # changed and unchanged
samplesH = np.zeros((len(terms),2 )) # changed and unchanged

knownparams = p.params["known"]
storage_home = os.getcwd() 

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 

pathname = path5

handles=[]
k = 0
Seeds = pickle.load(open(pathname+"Seeds.pickle","r"))
AllTrialsPD =[]
AllTrialsH=[]
indOrdTrialsPD=[]
indOrdTrialsH=[]
for trial,seed in enumerate(Seeds):
	for i,x in enumerate(terms):
		whichT = terms[i]
		print whichT	
		print x
		filename = "RHBM_PD_Self_Shuff_"+str(seed)+"_"+whichT+".pickle"
		if os.path.exists(pathname+filename)!=True:
			print pathname+filename+" not found !!!"
			continue
		RHBM_PD = pickle.load(open(pathname+filename,"r"))
		filename1 = "RHBM_H_Self_Shuff_"+str(seed)+"_"+whichT+".pickle"
		if os.path.exists(pathname1+filename1)!=True:
			print pathname1+filename1+" not found !!!"
			continue
		RHBM_H = pickle.load(open(pathname1+filename1,"r"))
		
		stillPD = len(RHBM_PD["stillPD"])
		PDchanged = len(RHBM_PD["changed"])
		stillH = len(RHBM_H["stillPD"])
		Hchanged = len(RHBM_H["changed"])
		print "stillPD",stillPD
		print "PDchanged",PDchanged
		print "stillH",stillH
		print "Hchanged",Hchanged

		samplesPD[i][0] = float(stillPD)/(stillPD+PDchanged)
		samplesPD[i][1] = float(PDchanged)/(stillPD+PDchanged)

		samplesH[i][0] = float(stillH)/(stillH+Hchanged)
		samplesH[i][1] = float(Hchanged)/(stillH+Hchanged)

		 
	AllTrialsPD.append(samplesPD)
	AllTrialsH.append(samplesH)

	# Sort the values tuple wise, that is highest amount of PD and healthy models unchanged will be replaced by its mean first. 
	tupleForm = [(x1,x2)   for x1,x2 in zip(samplesPD[:,0],samplesH[:,0])]
	tupleSorted = sorted(tupleForm,reverse=True)

	# Also sort PD separately and H separately
	PDSorted = sorted(samplesPD[:,0],reverse=True)
	HSorted = sorted(samplesH[:,0],reverse=True)
	indicesOrdPD = [ ]
	cnt =0
	for x in PDSorted:
		ind = np.where(samplesPD[:,0]==x)[0]
		if len(ind)>1 and cnt == 0:
			curr_ind=ind
			indicesOrdPD.append(ind[cnt])
			cnt+=1
		elif len(ind) >1 and cnt>0:
			if  len(np.where((curr_ind == ind)==False)[0]) == 0:
				indicesOrdPD.append(ind[cnt])
				cnt+=1
				if cnt == len(ind):
					cnt = 0
				else:
					cnt+=1
		else:
			indicesOrdPD.append(ind[0])
	
	indOrdTrialsPD.append(indicesOrdPD)
				
	indicesOrdH = [ ]   
	cnt =0
	for x in HSorted:
		ind = np.where(samplesH[:,0]==x)[0]
		if len(ind)>1 and cnt == 0:
			curr_ind=ind
			indicesOrdH.append(ind[cnt])
			cnt+=1
		elif len(ind) >1 and cnt>0:
			if  len(np.where((curr_ind == ind)==False)[0]) == 0:
				indicesOrdH.append(ind[cnt])
				cnt+=1
				if cnt == len(ind):
					cnt = 0
				else:
					cnt+=1
		else:
			indicesOrdH.append(ind[0])

	indOrdTrialsH.append(indicesOrdH)
	# Also find the corresponding indices
	indicesOrder=[]
	for i,x in enumerate(samplesPD[:,0]): # assuming that each tuple has an unique value of x, if this isnt true, the algortithm needs to change
	#	indicesOrder.append([y[0]   for y in tupleSorted].index(x))
		y = samplesH[i,0]
		for j,z in enumerate(tupleSorted):
			if z[0] == x and z[1] == y:
				indicesOrder.append(j)
				break


forReplaceInc=dict()
forReplaceInc["samplesPD"] = samplesPD
forReplaceInc["samplesH"] = samplesH
forReplaceInc["tupleSorted"] = tupleSorted
forReplaceInc["indicesOrder"] = indicesOrder
forReplaceInc["termsOrder"] = np.array(terms)[indicesOrder]
forReplaceInc["indicesOrdPD"] = indOrdTrialsPD
forReplaceInc["indicesOrdH"] = indOrdTrialsH
forReplaceInc["termsOrdPD"] = np.array(terms)[indOrdTrialsPD]
forReplaceInc["termsOrdH"] = np.array(terms)[indOrdTrialsH]	

pickle.dump(forReplaceInc,open(pathname+"forReplaceIncShuff.pickle","w"))



