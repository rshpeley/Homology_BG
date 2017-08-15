import numpy as np
params = dict()
#kscale=0.8
#kscale=0.75
kscale=0.75
params["known"] = {
	'd1d1' : -0.69,
	'd1d2' : -1.15,
	'd2d1' : -0.32,
	'd2d2' : -2.9,

	'd2fsi' : -0.53*0.6,
	'd1fsi' : -1.1*0.6,
	'gpid1' : -2.8,
	'gpiti' : -0.78,
	'gpistn' : 0.26,
	'stnstn' : 0.001*50,
	'gpita' : 0.0,
	'gpigpi' : 0.0,
	'fsifsi':-0.0012,

}
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
	'tistn' : np.linspace(0.01,posMax,noEl1)*kscale,
	'tastn' : np.linspace(0.01,posMax,noEl1)*kscale,
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


 
