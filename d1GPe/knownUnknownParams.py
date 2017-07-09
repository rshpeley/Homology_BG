import numpy as np
params = dict()
#kscale=0.8
#kscale=0.75
kscale=0.75
params["known"] = {
	#'d1d1' : -0.7*1.5,
	'd1d1' : -0.69,
	#'d1d2' : -0.73*1.5,
	'd1d2' : -1.15,
	#'d2d1' : -0.16*1.5,
	'd2d1' : -0.32,
	#'d2d2' : -0.97*1.5,
	'd2d2' : -2.9,

	#'d2fsi' : -0.22*1.5,
	'd2fsi' : -0.53*0.6,
	#'d1fsi' : -0.324*1.5,
	'd1fsi' : -1.1*0.6,
	#'gpid1' : -2.83*0.5,
	#'gpid1' : -21.,
	#'gpid1' : (-13.-2)/2.,
	'gpid1' : -2.8,
	#'gpid1' : -3.5,
#	'gpiti' : -0.17*2.5,
	#'gpiti' : -1.5,
	#'gpiti' : -0.34,
	#'gpiti' : -8.9,
	#'gpiti' : -3.9,
	#'gpiti' : -3.5,
	'gpiti' : -0.78,
	#'gpiti' : -0.5,
	#'gpita' : -0.17,
	#'gpistn' : 0.51,
	#'gpistn' : 0.8,
	#'gpistn' : 4.25,
	#'gpistn' : 3.8,
	
	#'gpistn' : 0.0,
	'gpistn' : 0.26,
	'stnstn' : 0.001*50,
	'gpita' : 0.0,
	'gpigpi' : 0.0,
	'fsifsi':-0.0012,
	'stnta':-0.0
#	'stnta' :-0.1,
#	'd1ti': -0.1,
#	'd2ti': -0.1

}
#params["iptau"] = 20.
params["iptau"] = 15.
negMax = -6.
posMax = 15.
noEl = 15
noEl1 = 20
# For lager time constants the weights had to be scaled down by 0.6
params["unknown"] = {
	'stnta' : np.linspace(-0.01,negMax,noEl)*kscale,
	'stnti' : np.linspace(-0.01,negMax,noEl)*kscale,
	#'tistn' : np.arange(0.01,1.25,0.25),
	'tistn' : np.linspace(0.01,negMax,noEl1)*kscale,
	'tastn' : np.linspace(0.01,negMax,noEl1)*kscale,
	'tata': np.linspace(-0.01,negMax,noEl)*kscale,
	'tati': np.linspace(-0.01,negMax,noEl)*kscale,
	'tita': np.linspace(-0.01,negMax,noEl)*kscale,
	'titi': np.linspace(-0.01,negMax,noEl)*kscale,
	'd1ta' : np.linspace(-0.01,negMax,noEl)*kscale,
	'd2ta' : np.linspace(-0.01,negMax,noEl)*kscale,
	'fsita' : np.linspace(-0.01,negMax,noEl)*kscale,
	'fsiti' : np.linspace(-0.01,negMax,noEl)*kscale,
	'tid2' : np.linspace(0,negMax,noEl)*kscale,
	'tad2' : np.linspace(0,negMax,noEl)*kscale,
	'tad1':np.linspace(0,negMax,noEl)*kscale,
	'tid1':np.linspace(0,negMax,noEl)*kscale,
	'd1ti':np.linspace(-0.01,negMax,noEl)*kscale,
	'd2ti':np.linspace(-0.01,negMax,noEl)*kscale,
	'jc1':np.linspace(1,posMax,noEl1)*kscale,
	'jc2':np.linspace(1,posMax,noEl1)*kscale,
	'jfsictx':np.linspace(1,posMax,noEl1)*kscale,
	'jstnctx':np.linspace(1,posMax,noEl1)*kscale

}
params["dt"] = 0.01
# For transfer function, since the sampling size is already in msecs, dt = 1, dt ~ nyquist frequency = 2*maximum frequency in signal, 
# maximum frequency that can be detected is 500 hz 
params["dtTF"] = 1


 
