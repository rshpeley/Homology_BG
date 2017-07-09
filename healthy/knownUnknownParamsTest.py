import numpy as np
params = dict()
#kscale=0.8
kscale=0.75
params["known"] = {
	'd1d1' : -0.7*1.5,
	'd1d2' : -0.73*1.5,
	'd2d1' : -0.16*1.5,
	'd2d2' : -0.97*1.5,

	'd2fsi' : -0.22*1.5,
	'd1fsi' : -0.324*1.5,
	#'gpid1' : -2.83*0.5,
	'gpid1' : -2.83*1.5,
#	'gpiti' : -0.17*2.5,
	#'gpiti' : -1.5,
	'gpiti' : -0.34,
	'gpita' : -0.17,
	#'gpistn' : 0.51,
	'gpistn' : 0.8,
	'stnstn' : 0.001*50,
	'gpita' : 0.0,
	'gpigpi' : 0.0,
	'fsifsi':-0.0012,
#	'stnta' :-0.1,
#	'd1ti': -0.1,
#	'd2ti': -0.1

}
# For lager time constants the weights had to be scaled down by 0.6
params["unknown"] = {
	'stnta' : np.arange(-0.01,-3.,-0.6)*kscale,
	'stnti' : np.arange(-0.01,-3.,-0.6)*kscale,
	#'tistn' : np.arange(0.01,1.25,0.25),
	'tistn' : np.arange(0.01,3.,0.6)*kscale,
	'tastn' : np.arange(0.01,3.,0.6)*kscale,
	'tata': np.arange(-0.01,-2.,-0.4)*kscale,
	'tati': np.arange(-0.01,-2.,-0.4)*kscale,
	'tita': np.arange(-0.01,-2.,-0.4)*kscale,
	'titi': np.arange(-0.01,-2.,-0.4)*kscale,
	'd1ta' : np.arange(-0.01,-2.,-0.4)*kscale,
	'd2ta' : np.arange(-0.01,-2.,-0.4)*kscale,
	'fsita' : np.arange(-0.01,-2.,-0.4)*kscale,
	'fsiti' : np.arange(-0.01,-2.,-0.4)*kscale,
	'tid2' : np.arange(0,-5,-1)*kscale,
	'tad2' : np.arange(0,-5,-1)*kscale,
	'd1ti':np.arange(-0.01,-2.,-0.4)*kscale,
	'd2ti':np.arange(-0.01,-2.,-0.4)*kscale,
	'jc1':np.arange(1,15,3)*kscale,
	'jc2':np.arange(1,15,3)*kscale,
	'jfsictx':np.arange(1,15,3)*kscale,
	'jstnctx':np.arange(1,20,4)*kscale

}
params["dt"] = 0.01
params["dtTF"] = 1 # Check ../knownUnknown for why 


params1 = dict()
kscale=1
params1["known"] = {
	'd1d1' : -0.7*1.,
	'd1d2' : -0.73*1.,
	'd2d1' : -0.16*1.,
	'd2d2' : -0.97*1.,

	'd2fsi' : -0.22*1.,
	'd1fsi' : -0.324*1.,
	#'gpid1' : -2.83*0.5,
	'gpid1' : -2.83,
#	'gpiti' : -0.17*2.5,
	#'gpiti' : -1.5,
	'gpiti' : -1.0,
	'gpita' : -0.17*kscale,
	#'gpistn' : 0.51,
	'gpistn' : 0.51,
	'stnstn' : 0.001*10,
	'gpita' : 0.0,
	'gpigpi' : 0.0,
	'fsifsi':-0.0012*kscale,
}

params1["unknown"] = {
	'stnta' : np.arange(-0.01,-0.65,-0.15),
	'stnti' : np.arange(-0.01,-0.65,-0.15),
	'tistn' : np.arange(0.01,1.25,0.25),
	'tastn' : np.arange(0.01,1.25,0.25),
	'tata': np.arange(-0.01,-0.65,-0.15),
	'tati': np.arange(-0.01,-0.65,-0.15),
	'tita': np.arange(-0.01,-0.65,-0.15),
	'titi': np.arange(-0.01,-0.65,-0.15),
	'd1ta' : np.arange(-0.01,-0.65,-0.15),
	'd2ta' : np.arange(-0.01,-0.65,-0.15),
	'fsita' : np.arange(-0.01,-0.65,-0.15),
	'fsiti' : np.arange(-0.01,-0.65,-0.15),
	'tid2' : np.arange(0,-5,-1),
	'tad2' : np.arange(0,-5,-1),
	'd1ti':np.arange(-0.01,-0.65,-0.15),
	'd2ti':np.arange(-0.01,-0.65,-0.15),
	'jc1':np.arange(1,6.6,1.3),
	'jc2':np.arange(1,6.6,1.3),
	'jfsictx':np.arange(1,6.6,1.3),
	'jstnctx':np.arange(1,6.6,1.3)
}

params1["dt"] = 0.01

'''
params["unknown"] = {
	'stnta' : np.arange(-0.01,-2,-0.18),
	'stnti' : np.arange(-0.01,-2,-0.18),
	'tistn' : np.arange(0.01,3,0.25),
	'tastn' : np.arange(0.01,3,0.25),
	'tata': np.arange(-0.01,-2,-0.18),
	'tati': np.arange(-0.01,-2,-0.18),
	'tita': np.arange(-0.01,-2,-0.18),
	'titi': np.arange(-0.01,-2,-0.18),
	'd1ta' : np.arange(-0.01,-2,-0.18),
	'd2ta' : np.arange(-0.01,-2,-0.18),
	'fsita' : np.arange(-0.01,-3,-0.25),
	'fsiti' : np.arange(-0.01,-3,-0.25),
	'tid2' : np.arange(0,-3,-0.25),
	'tad2' : np.arange(0,-3,-0.25),
	'd1ti':np.arange(-0.01,-2,-0.18),
	'd2ti':np.arange(-0.01,-2,-0.18),
	'jc1':np.arange(0.5,8,0.65),
	'jc2':np.arange(0.5,8,0.65),
	'jfsictx':np.arange(0.5,8,0.65),
	'jstnctx':np.arange(0.5,8,0.65)

}
''' 
