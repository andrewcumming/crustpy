from scipy.optimize import brentq
import numpy
from scipy import interpolate
from math import exp, pi, sqrt, log

class Eos:
	def __init__(self,T=1e8,P=1e27,A=56,Z=26,Yn=0.0,Qimp=1):
		self.P=P
		self.T=T
		self.A=A
		self.Z=Z
		self.Yn=Yn
		self.Ye=(1.0-Yn)*(1.0*Z/A)
		self.Acell=A/(1.0-Yn)
		self.Yi=1.0/self.Acell
		self.Qimp=Qimp
		self.calculate_rho()
		
	def __str__(self):
		return "P=%g rho=%g T=%g A=%g Z=%g Qimp=%g Yn=%g EFermi=%g ne=%g Ye=%g Acell=%g Tp=%g" % (self.P,self.rho,self.T,self.A,self.Z,self.Qimp,self.Yn,self.EFermi,self.ne,self.Ye,self.Acell,self.TP)

	def calculate_rho(self):
		self.rho = brentq(self.pressure,2e7,2e14,xtol=1e-6)
		self.ne = self.rho*self.Ye/1.67e-24
		self.EFermi = 6.09e-11*self.ne**(1.0/3.0)     # T=0 EF in MeV
		self.x=self.EFermi/0.511
		# In the next line, we are using self.A as the mass of the nucleus (no entrainment)
		self.TP=7.76e3*self.Z*sqrt(self.Yi*self.rho/self.A)

	def pressure(self,rho):
		# used in the root-find to find the density
		self.rho=rho
		Pe=1.231e15*(self.rho*self.Ye)**(4.0/3.0)
		# Neutron pressure from Mackie & Baym 
		k=0.207*(1e-12*self.rho*self.Yn)**(1.0/3.0)
		EFn=1.730*k+25.05*k*k-30.47*k*k*k+17.42*k*k*k*k
		Pn=0.4*EFn*1.6e-6*self.rho*self.Yn/1.67e-24
		return Pe + Pn - self.P
		
	def update_T(self,T):
		self.T=T

	def Kcond(self):
		# Thermal conductivity
		# Phonon scattering frequency
		lambda_ep=1.0		
		TU=8.7*self.rho**0.5*(self.Ye/0.05)*(self.Z/30.0)**(1.0/3.0)
		TD=3.5e3*self.Ye*self.rho**0.5
		fep=1.247e10*self.T*lambda_ep*exp(-TU/self.T)/(1.0+(TD/(3.5*self.T))**2)**0.5
		# Impurity scattering frequency
		lambda_eQ=0.5*log(1.0+0.4*137.0*pi)*(1.0+2.5/(137.0*pi))-0.5
		feQ=1.75e16*self.x*self.Qimp*lambda_eQ/self.Z
		fc = fep+feQ
		return 4.116e19*self.T*self.rho*self.Ye/(self.x*fc)
				
	def CV(self):
		# Heat capacity has contributions from electrons, lattice, and neutrons
		return self.CV_electrons()+self.CV_ions()+self.CV_neutrons()
		
	def CV_electrons(self):
		# Electron contribution to the heat capacity
		return 7.12e-3*pi**2*self.Ye*self.T/self.EFermi

	def CV_neutrons(self):
		# Neutron contribution to the heat capacity
		return 0.0

	def CV_ions(self):
		# Ion contribution to the heat capacity	
		# from equation (5) of Chabrier 1993
		eta=self.TP/self.T
		# x is alpha eta,  y is gamma eta
		x=0.399*eta
		y=0.899*eta
		ex=exp(x)
		ey=exp(y)
		dd1=pi**4/(5*x**3) - 3.0*(6.0+x*(6.0+x*(3.0+x)))/(ex*x**3)
		dd2=1.0-0.375*x+0.05*x**2
		dd=min(dd1,dd2)
		return 8.26e7*self.Yi*(8.0*dd-6*x/(ex-1.0)+(y**2*ey/(ey-1.0)**2))
		
		
		
if __name__ == '__main__':
	myeos = Eos()
	print myeos
	