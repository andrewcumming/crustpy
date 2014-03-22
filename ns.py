from scipy.optimize import brentq

def grav(M,R):
	# returns g14 and 1+z
	g = (6.67e-8*2e33/(1e14*1e10)) * M / R**2
	f = (2.0*1e14*g*1e5*R/9e20)
	zz = 1.0
	if f<1.0:
		zz = 1.0/(1.0-f)**0.5
		g*=zz
	else:
		g=1e10
	return g,zz

def grav_eqn(R,M,g):
	return grav(M,R) - g

def R(M,g):
	return brentq(grav_eqn,1.0,20.0,xtol=1e-6,args=(M,g,))
