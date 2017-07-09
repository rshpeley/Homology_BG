import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
#mpl.use('Agg')
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
import sys
import paramsearchGA_DopDep_nonlinear as psGA
import knownUnknownParams as p


pathname = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/output/"
#PDFeatures = pickle.load(open(pathname+"PDFeaturesnew.pickle","r"))
ipNo1 = 5
PDFeatures = pickle.load(open(pathname+"PDFeaturesnew_"+str(ipNo1)+".pickle","r"))



clustersPD = pickle.load(open(pathname+"clustersPDMore"+str(ipNo1)+".pickle","r"))

labelPD = pickle.load(open(pathname+"LabeledClustersPD"+str(ipNo1)+".pickle","r"))


gmmPD = pickle.load(open(pathname+"gmmPDwo"+str(ipNo1)+".pickle","r"))


knownparams = p.params["known"]

ipctx1=dict()
ipctx1["ip"] = np.zeros((1,2001))


RatesD1 = dict()
RatesD2 = dict()


def calcFI(params):
	for i,clus1 in enumerate(clustersPD[::1]):
		if len(clus1) > 1:
			# Extract the means from GMMs
			print i
			means = gmmPD[i]["means"]
			meanD1 =[]
			meanD2 = []
			ipamp = params["ipAmp"]
			targLabel = labelPD[i]["label"]
			corrInds = []
			allLabels=[]
			for j,ind in enumerate(means):
			#if gmmPD[i]["labels"][j]["label"] == targLabel:
				# Construct A and B out of mean
				d1ta = ind[0]
				d2ta = ind[1]
				fsita = ind[2]
				fsiti = ind[3]
				tata = ind[4]
				tati = ind[5]
				tastn = ind[6]
				tita = ind[7]
				titi = ind[8]
				tistn = ind[9]
				stnta = ind[10]
				stnti = ind[11]			
				tid2 = ind[12]
				tad2 = ind[13]	
				d1ti = ind[14]
				d2ti = ind[15]
				jc1 = ind[16]
				jc2 = ind[17]
				#jc2 = jc1 # Prevent trivial solutions
				jfsictx = ind[18]
				jstnctx = ind[19]

				A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,d1ti,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,d2ti,0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,tad2,0.,tata,tati,tastn,0],[0.,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,stnta,stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
				#print A
				B = np.matrix([jc1,jc2,jfsictx,0,0,jstnctx,0])
				Flags = []
				Flags.append("Trans")

				Rates = psGA.calcRates(Flags,1,A,B,False,ipctx1,ipamp=ipamp,iptau=14.)

				meanD1.append(np.mean(Rates['d1']))
				meanD2.append(np.mean(Rates['d2']))
				corrInds.append(j)
				allLabels.append(gmmPD[i]["labels"][j]["label"])
			RatesD1[i] =dict()
			#RatesD1[i]["label"] = labelPD[i]["label"]
			RatesD1[i]["label"] = allLabels
			RatesD1[i]["F_I"] = meanD1
			RatesD1[i]["corrInds"]=corrInds
			RatesD2[i] =dict()
			#RatesD2[i]["label"] = labelPD[i]["label"]
			RatesD2[i]["label"] = allLabels
			RatesD2[i]["F_I"] = meanD2
			RatesD2[i]["corrInds"]=corrInds
	pickle.dump(RatesD1,open(pathname+"FI_clusterWise_D1_"+str(ipamp)+".pickle","w"))
	pickle.dump(RatesD2,open(pathname+"FI_clusterWise_D2_"+str(ipamp)+".pickle","w"))
