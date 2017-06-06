#cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True
from scipy.optimize import brentq
#from numpy import pi
from libc.math cimport exp, log, sqrt

def init_eos(double T=1e8,double P=1e27,double A=56,double Z=26,double Yn=0.0,double Qimp=1):
	cdef double Acell = A/(1.0-Yn)
	eos = dict(T=T,P=P,A=A,Z=Z,Yn=Yn,Qimp=Qimp,Ye=(1.0-Yn)*(1.0*Z/A),\
		Acell=Acell,Yi=1.0/Acell)
	eos = calculate_rho(eos)
	eos = update_T(eos,T)
	return eos

def calculate_rho(eos):
	cdef double rho, ne, EFermi, x, TP, Yi, A, Ye, Z, Yn, P
	A = eos['A']
	Z = eos['Z']
	Yi = eos['Yi']	
	Ye = eos['Ye']	
	Yn = eos['Yn']	
	P = eos['P']	
	rho = brentq(pressure,2e7,2e14,xtol=1e-6,args=(Yn,Ye,P,))
	ne = rho*Ye/1.67e-24
	EFermi = 6.09e-11*ne**(1.0/3.0)  # T=0 EF in MeV
	x = EFermi/0.511
	# In the next line, we are using self.A as the mass of the nucleus (no entrainment)
	TP=7.76e3*Z*sqrt(Yi*rho/A)
	eos['rho'] = rho
	eos['ne'] = ne
	eos['x'] = x
	eos['TP'] = TP
	eos['EFermi'] = EFermi	
	return eos
 	
cdef double pressure(double rho,double Yn, double Ye, double P):
	#used in the root-find to find the density
	cdef double Pe, k, EFn, Pn
	Pe=1.231e15*(rho*Ye)**(4.0/3.0)
	# Neutron pressure from Mackie & Baym 
	k=0.207*(1e-12*rho*Yn)**(1.0/3.0)
	EFn=1.730*k+25.05*k*k-30.47*k*k*k+17.42*k*k*k*k
	Pn=0.4*EFn*1.6e-6*rho*Yn/1.67e-24
	return Pe + Pn - P
		
def update_T(eos,double T):
	cdef double TP = eos['TP']
	cdef double Ye = eos['Ye']
	cdef double Yi = eos['Yi']
	cdef double x = eos['x']
	cdef double Z = eos['Z']
	cdef double EFermi = eos['EFermi']
	cdef double Qimp = eos['Qimp']
	cdef double rho = eos['rho']
	eos['T']=T
	eos['Kcond'] = calculate_Kcond(T,Ye,Z,x,Qimp,rho)
	eos['CV_electrons'] = CV_electrons(T,Ye,EFermi)
	eos['CV_ions'] = CV_ions(T,TP,Yi)
	eos['CV_neutrons'] = CV_neutrons(T)
	eos['CV'] = eos['CV_electrons'] + eos['CV_ions'] + eos['CV_neutrons']
	#eos['CV'] = calculate_CV(T,TP,Ye,Yi,EFermi)
	return eos

cdef double calculate_Kcond(double T,double Ye,double Z,double x,double Qimp,double rho):
	# Thermal conductivity
	cdef double lambda_ep
	cdef double TU, TD, fep, lambda_eQ, feQ, fc
	cdef double mypi = 3.141592654
	# Phonon scattering frequency
	lambda_ep=1.0		
	TU=8.7*sqrt(rho)*(Ye/0.05)*(Z/30.0)**(1.0/3.0)
	TD=3.5e3*Ye*sqrt(rho)
	fep=1.247e10*T*lambda_ep*exp(-TU/T)/sqrt(1.0+(TD/(3.5*T))**2)
	# Impurity scattering frequency
	lambda_eQ=0.5*log(1.0+0.4*137.0*mypi)*(1.0+2.5/(137.0*mypi))-0.5
	feQ=1.75e16*x*Qimp*lambda_eQ/Z
	fc = fep+feQ
	return 4.116e19*T*rho*Ye/(x*fc)
				
cdef double calculate_CV(double T,double TP,double Ye,double Yi,double EFermi):
	# Heat capacity has contributions from electrons, lattice, and neutrons
	return CV_electrons(T,Ye,EFermi) + CV_ions(T,TP,Yi) + CV_neutrons(T)
	
cdef double CV_electrons(double T, double Ye, double EFermi):
	# Electron contribution to the heat capacity
	cdef double mypi = 3.141592654
	return 7.12e-3*mypi**2*Ye*T/EFermi
	
cdef double CV_neutrons(double T):
	# Neutron contribution to the heat capacity
	return 0.0

cdef double CV_ions(double T, double TP, double Yi):
	# Ion contribution to the heat capacity	
	# from equation (5) of Chabrier 1993
	cdef double eta, x, y, ex, ey, dd1, dd2, dd, fac
	cdef double mypi = 3.141592654
	eta=TP/T
	# x is alpha eta,  y is gamma eta
	x=0.399*eta
	y=0.899*eta
	ex=exp(x)
	ey=exp(y)
	dd1=mypi**4/(5*x**3) - 3.0*(6.0+x*(6.0+x*(3.0+x)))/(ex*x**3)
	dd2=1.0-0.375*x+0.05*x**2
	dd=min(dd1,dd2)
	fac=8.0*dd-6*x/(ex-1.0)+(y**2*ey/(ey-1.0)**2)
	return 8.26e7*Yi*fac
