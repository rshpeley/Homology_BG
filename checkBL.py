import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')


import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import scipy.fftpack as spfft
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *
#import paramsearchGA_DopDep as psGA
import paramsearchGA_DopDep_nonlinear as psGA
import knownUnknownParams as p

AllAs=[]
AllASTPs=[]
AllBs=[]
storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams'
from scipy.optimize import fsolve
ctx = 5.0
terms1 = ["d1","d2","fsi","ta","ti","stn","gpi","ipctx"]



def run(params):
	pathname = storage_home+'/output/'
	pathname1 = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams'+'/output/'
	print "in classification"
	knownparams = p.params["known"] 
	leak = -0.05
	Supressed=[]
	Sup_Ind=[]
	k = 0
	randSamps = pickle.load(open(pathname+"randSampsTest.pickle","r"))
	Allcombs = pickle.load(open(pathname+"AllcombsTest.pickle","r"))


	ipAmp = params["ipAmp"]

	print ipAmp
	baseLineFR = []
	for i,j in enumerate(randSamps):
		print i
		ind = Allcombs[j]
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
		jfsictx = ind[18]
		jstnctx = ind[19]

		#A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,knownparams['d1ti'],0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,knownparams['d2ti'],0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,tad2,0.,tata,tati,tastn,0],[0.,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,knownparams['stnta'],stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
		A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,d1ti,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,d2ti,0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,tad2,0.,tata,tati,tastn,0],[0.,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,stnta,stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])

		B = np.matrix([jc1,jc2,jfsictx,0,0,jstnctx,0])
		ipctx1=dict()
		#Calculate Rates for SWA and lesion(dopamine depletion)
		Flags = []
		Flags.append("Trans")
		#Flags.append("DopDep")
		
		ipctx1["ip"] = np.zeros((1,2001))
		TransRates = psGA.calcRates(Flags,1.,A,B,False,ipctx1,ipAmp,iptau=p.params["iptau"])
		temp = dict()
		for x in terms1:
			temp[x] = np.mean(TransRates[x])
		baseLineFR.append(temp) 

	pickle.dump(baseLineFR,open(pathname1+"baseLineFR_"+str(ipAmp)+".pickle","w"))
