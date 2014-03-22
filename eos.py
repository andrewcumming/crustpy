from scipy.optimize import brentq

class Eos:
	def __init__(self,T=1e8,P=1e27,A=56,Z=26,Yn=0.0,Qimp=1):
		self.T=T
		self.P=P
		self.A=A
		self.Z=Z
		self.Yn=Yn
		self.Ye=(1.0-Yn)*(1.0*Z/A)
		self.Qimp=Qimp
		self.calculate_rho()
		
	def __str__(self):
		return "P=%g rho=%g T=%g A=%g Z=%g Qimp=%g Yn=%g EFermi=%g ne=%g Ye=%g" % (self.P,self.rho,self.T,self.A,self.Z,self.Qimp,self.Yn,self.EFermi,self.ne,self.Ye)

	def calculate_rho(self):
		self.rho = brentq(self.pressure,2e7,2e14,xtol=1e-6)
		self.ne = self.rho*self.Ye/1.67e-24
		self.EFermi = 6.09e-11*self.ne**(1.0/3.0)     # in MeV

	def pressure(self,rho):
		self.rho=rho
		Pe=1.231e15*(self.rho*self.Ye)**(4.0/3.0)
		Pn=0.0
		if self.Yn > 0.0:
			k=0.207*(1e-12*self.rho*self.Yn)**(1.0/3.0)
			EFn=1.730*k+25.05*k*k-30.47*k*k*k+17.42*k*k*k*k
			Pn=0.4*EFn*1.6e-6*self.rho*self.Yn/1.67e-24
		return Pe + Pn - self.P
		
	def update_T(self,T):
		self.T=T

	def Kcond(self):
		# Thermal conductivity
		return 1e21

	def CV(self):
		# Heat capacity
		return 1e21
		
		
if __name__ == '__main__':
	myeos = Eos()
	print myeos
	