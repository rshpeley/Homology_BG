import pickle
import numpy as np
import pylab as pl

baseLine = pickle.load(open("/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/output/baseLineFR_0.pickle","r"))
Rates=dict()
nuclei = [x for x in baseLine[0].keys()]
for x in nuclei:
    temp = []                                     
    for el in baseLine:
        temp.append(el[x])
    Rates[x]=temp

lowbins = np.arange(0,50,0.5)
Hists=dict()
pl.figure()
for i,x in enumerate(nuclei):
    Hists[x] = np.histogram(Rates[x],bins=lowbins)
    pl.plot(Hists[x][1][:-1],Hists[x][0],'-',label = x)

pl.legend()
pl.show()	

