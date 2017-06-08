#cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True

from scipy.optimize import brentq

cdef extern from "math.h":
	double exp(double m)
	double log(double m)
	double sqrt(double m)
	double MPI "M_PI"
	        
cdef extern:
	void condegin_(double *temp,double *densi,double *B,double *Zion,double *CMI,
			double *CMI1,double *Zimp, double *RSIGMA,double *RTSIGMA,
			double *RHSIGMA,double *RKAPPA,double *RTKAPPA,double *RHKAPPA)

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
	cdef double Yn = eos['Yn']
	cdef double A = eos['A']
	cdef double x = eos['x']
	cdef double Z = eos['Z']
	cdef double EFermi = eos['EFermi']
	cdef double Qimp = eos['Qimp']
	cdef double rho = eos['rho']
	eos['T']=T
	# simpler version of conductivity that has lambda_ep set to 1
	# (gives ~factor of 2 difference in derived Qimp for MXB 1659)
	#eos['Kcond'] = calculate_Kcond(T,Ye,Z,x,Qimp,rho)
	eos['Kcond'] = potek_cond(rho,T,A,Z,Yn,Qimp)
	eos['CV_electrons'] = CV_electrons(T,Ye,EFermi)
	eos['CV_ions'] = CV_ions(T,TP,Yi)
	eos['CV_neutrons'] = CV_neutrons(rho, T, Yn)
	eos['CV'] = eos['CV_electrons'] + eos['CV_ions'] + eos['CV_neutrons']
	#eos['CV'] = calculate_CV(rho, T,TP,Ye,Yi,EFermi, Yn)
	return eos

cdef double potek_cond(double rho, double T, double A, double Z, double Yn, double Qimp):
	# returns the thermal conductivity in cgs from Potekhin's fortran code
	cdef double s1,s2,s3,k1,k2,k3
	cdef double Zimp = sqrt(Qimp)
	cdef double Acell = A / (1.0-Yn)
	cdef double Bfield = 0.0
	cdef double temp = T*1e-6/5930.0
	cdef double rr = rho/(Acell*15819.4*1822.9)
	condegin_(&temp,&rr,&Bfield,&Z,&A,&Acell,&Zimp, &s1,&s2,&s3,&k1,&k2,&k3)
	return k1*2.778e15
	
cdef double calculate_Kcond(double T,double Ye,double Z,double x,double Qimp,double rho):
	# Thermal conductivity
	cdef double lambda_ep
	cdef double TU, TD, fep, lambda_eQ, feQ, fc
	# Phonon scattering frequency
	lambda_ep=1.0		
	TU=8.7*sqrt(rho)*(Ye/0.05)*(Z/30.0)**(1.0/3.0)
	TD=3.5e3*Ye*sqrt(rho)
	fep=1.247e10*T*lambda_ep*exp(-TU/T)/sqrt(1.0+(TD/(3.5*T))**2)
	# Impurity scattering frequency
	lambda_eQ=2.09 # = 0.5*log(1.0+0.4*137.0*MPI)*(1.0+2.5/(137.0*MPI))-0.5
	feQ=1.75e16*x*Qimp*lambda_eQ/Z
	fc = fep+feQ
	return 4.116e19*T*rho*Ye/(x*fc)
				
cdef double calculate_CV(double rho, double T,double TP,double Ye,double Yi,double EFermi, double Yn):
	# Heat capacity has contributions from electrons, lattice, and neutrons
	return CV_electrons(T,Ye,EFermi) + CV_ions(T,TP,Yi) + CV_neutrons(rho, T, Yn)
	
cdef double CV_electrons(double T, double Ye, double EFermi):
	# Electron contribution to the heat capacity
	return 7.12e-3*MPI**2*Ye*T/EFermi
	
cdef double CV_ions(double T, double TP, double Yi):
	# Ion contribution to the heat capacity	
	# from equation (5) of Chabrier 1993
	cdef double eta, x, y, ex, ey, dd1, dd2, dd, fac
	eta=TP/T
	# x is alpha eta,  y is gamma eta
	x=0.399*eta
	y=0.899*eta
	ex=exp(x)
	ey=exp(y)
	dd1=MPI**4/(5*x**3) - 3.0*(6.0+x*(6.0+x*(3.0+x)))/(ex*x**3)
	dd2=1.0-0.375*x+0.05*x**2
	dd=min(dd1,dd2)
	fac=8.0*dd-6*x/(ex-1.0)+(y**2*ey/(ey-1.0)**2)
	return 8.26e7*Yi*fac

cdef double CV_neutrons(double rho, double T, double Yn):
	# Neutron contribution to the heat capacity
	cdef double k, EFn, cvneut
	cdef double R00, t, u

	if Yn == 0.0:
		return 0.0
	# Mackie and Baym
	k = 0.207*(1e-12*rho*Yn)**(1.0/3.0)
	EFn = 1.730*k + 25.05*k*k - 30.47*k*k*k + 17.42*k*k*k*k # in MeV
	cvneut = 0.5*3.1415*3.1415*8.3144e7*Yn * 1.38e-16*T/(EFn*1.6e-6)

	# gap suppression
	t = T / TC(rho, T, Yn)	
	if t>1.0:
		R00 = 1.0
	else:
		u = sqrt(1.0-t)*(1.456 - 0.157/sqrt(t)+1.764/t)
		R00 = (0.4186 + sqrt(1.007*1.007 + (0.5010*u)**2.0))**2.5 * exp(1.456 - sqrt(1.456*1.456+u*u))
	if R00 < 1e-8:
		R00 = 0.0
	cvneut *= R00

	return cvneut

cdef double TC(double rho, double T, double Yn):
	# calculates the critical temp for neutron superfluidity in the crust
	cdef double k, Tcrit, delk, a, b, t
	cdef int i1, i2, i
	cdef double k0[17]
	cdef double d0[17]
	cdef double d2[17]
	
	# neutron k_F in fm^-1
	k = 0.261*(1e-12*rho*Yn)**(1.0/3.0)
	
	# gap from SFB03
	k0 = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.175,1.25, 1.3, 1.35, 1.4, 1.45)
	d0 = (0.0, 0.09, 0.210, 0.360, 0.500, 0.610, 0.720, 0.790, 0.780,0.700, 0.580, 0.450, 0.280, 0.190, 0.100, 0.030, 0.0)
	d2 = (2.915e1, -4.297, 6.040, -1.863, -4.59, 2.221, -4.296,-9.037, -7.555, -2.741, -5.480, -1.344e1, 1.656e1, -6.667,1.010e1, 1.426e1, 2.887e1)
	if k < k0[0] or k > k0[16]:
		return 0.0
	i1=0
	i2=16
	while (i2-i1)>1:
		i = (i2+i1)/2
		if k0[i]>k:
			i2 = i
		else:
			i1 = i
	delk = k0[i2]-k0[i1]
	a = (k0[i2]-k)/delk
	b = (k-k0[i1])/delk
	t = a*d0[i1] + b*d0[i2] + (((a**3.0)-a)*d2[i1] + ((b**3.0)-b)*d2[i2]) * (delk*delk)/6
	
	Tcrit = (t/1.76) * 1.1604e10
	return Tcrit
