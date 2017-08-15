import pickle
import knownUnknownParamsStrRed as p
import numpy as np

uniqCombs = pickle.load(open("output/uniqueCombs.pickle","r"))
knownparams = p.params["known"]

AllAs = []
AllBs = []

for ind in uniqCombs:
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

	A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,d1ti,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,d2ti,0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,tad2,0.,tata,tati,tastn,0],[0.,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,stnta,stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
	B = np.matrix([jc1,jc2,jfsictx,0,0,jstnctx,0])

	AllAs.append(A)
	AllBs.append(B)

pickle.dump(AllAs,open("output/AllAs.pickle","w"))
pickle.dump(AllBs,open("output/AllAs.pickle","w"))

