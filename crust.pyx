#cython: boundscheck=False, wraparound=True, nonecheck=False, cdivision=True
import ns
import eos
import numpy
from scipy.integrate import odeint
from scipy import interpolate
#from numpy import log, exp, log10
from libc.math cimport exp, log, sqrt, log10

class Crust:
	def __init__(self,mass=1.4,radius=12.0,ngrid=30,Tc=1e8,Qimp=1,shallow_heating=False,shallow_y=1e13,shallow_Q=1.0,cooling_bc=False):
		self.mass = mass
		self.radius = radius
		self.g,self.zz = ns.grav(mass,radius)
		self.g *= 1e14
		self.ngrid=ngrid
		self.Tc = Tc
		self.Qimp = Qimp
		self.shallow_heating = shallow_heating
		self.shallow_Q = shallow_Q
		self.shallow_y = shallow_y
		self.cooling_bc = cooling_bc
		P1 = 1e12*2.28e14
		P2 = 6.5e32
		self.grid = []
		self.P=[]
		self.dx = log(P2/P1)/(ngrid-1)
		for i in range(0,ngrid):
			this_P = exp(log(P1) + self.dx*i)
			(A,Z,Yn) = self.composition(this_P)
			grid_point = eos.init_eos(P=this_P,Yn=Yn,A=A,Z=Z,Qimp=Qimp,T=Tc)
			self.grid.append(grid_point)
			self.P.append(this_P)
		self.P=numpy.array(self.P)
		# Read in the envelope Tb-Teff relation
		F = numpy.array([float(line.split()[3]) for line in open('TbTeff','r')])
		T = numpy.array([float(line.split()[2]) for line in open('TbTeff','r')])
		self.TbTeff = interpolate.interp1d(T,F)
		#self.pool = Pool(4)
								
	def __str__(self):
		return "M=%g R=%g g=%g ngrid=%d Tc=%g Qimp=%g" % (self.mass,self.radius,self.g/1e14,self.ngrid,self.Tc,self.Qimp)

	def temperature_profile(self):
		return numpy.array([item['rho'] for item in self.grid]),numpy.array([item['T'] for item in self.grid])

	def set_top_temperature(self,T):
		self.grid[0] = eos.update_T(self.grid[0],T)

	def set_core_temperature(self,T):
		self.grid[-1] = eos.update_T(self.grid[-1],T)

	def composition(self,double P):
		cdef double Yn, Acell, A, Z, Pmax
		# composition as a function of pressure from HZ 1990
		Acella=numpy.array([56.0,56.0,56.0,56.0,56.0,56.0,56.0,56.0,112.0,112.0,112.0,112.0,112.0,224.0,224.0,224.0,224.0,448.0,448.0])
		Aa=numpy.array([56.0,56.0,56.0,56.0,56.0,52.0,46.0,40.0,68.0,62.0,56.0,50.0,44.0,66.0,60.0,54.0,48.0,96.0,88.0]);
		Za=numpy.array([26.0,24.0,22.0,20.0,18.0,16.0,14.0,12.0,20.0,18.0,16.0,14.0,12.0,18.0,16.0,14.0,12.0,24.0,22.0]);
		# Maximum pressure at which that composition is present
		Pmaxa = numpy.array([7.235e26,9.569e27,1.152e29,4.747e29,1.361e30,1.980e30,2.253e30,2.637e30,2.771e30,3.216e30,3.825e30,4.699e30,6.043e30,7.233e30,9.238e30,1.228e31,1.602e31,1.613e31,1e33])
		# Find the composition in the first entry in the table where P>Pmax.
		(Acell,A,Z) = [(Acell,A,Z) for (Acell,A,Z,Pmax) in zip(Acella,Aa,Za,Pmaxa) if Pmax>P][0]
		Yn = (Acell-A)/Acell
		return (A,Z,Yn)

	def evolve(self,time=10000.0,mdot=0.1,n=100):
		print("Evolving in time for %g days at mdot=%g" % (time,mdot))
		inic = numpy.array([item['T'] for item in self.grid])
		times = 10.0**(numpy.arange(n)*log10(time/self.zz)/(n-1))
		self.eps = numpy.array([mdot*self.crust_heating(P) for P in self.P])
		result,info = odeint(self.derivs,inic,times*3600.0*24.0,rtol=1e-6,atol=1e-6,args=(mdot,),full_output=True)
		for i,T in enumerate(result[-1,:]): 
			self.grid[i] = eos.update_T(self.grid[i],T)
		Teff = [(self.envelope_flux(T)/5.67e-5)**0.25*1.38e-16/(1.6e-12*self.zz) for T in result[:,0]]
		# returns the redshifted lightcurve
		return times*self.zz, Teff
		
#	def update_T_j(self,args):
#		#print(args)
#		j, T, grid = args
#		#print("hello from update_T:",j,T)
#		grid.update_T(T)
#		return True
			
	def derivs(self,T,double time,double mdot):
		# calculate fluxes
		#result = self.pool.map(self.update_T_j, zip(range(self.ngrid),T,self.grid))
		#list(map(self.update_T_j,zip(range(self.ngrid),T)))
		self.grid = [eos.update_T(grid_j,T_j) for grid_j,T_j in zip(self.grid,T)]
	#	for j in range(0,self.ngrid):
	#		self.grid[j] = eos.update_T(self.grid[j],T[j]) 
		vals = numpy.array([[item['rho']*self.g*item['Kcond']/item['P'],item['CV']] for item in self.grid])
		K=vals[:,0]
		cv=vals[:,1]
		# F[j] is the flux at j-1/2
		F=numpy.array([0.0]*self.ngrid)
		F[1:] = 0.5*(K[:-1]+K[1:])*(T[1:]-T[:-1])
		F/=self.dx
		if mdot>0.0 and not self.cooling_bc:
			F[0]=F[1]
		else:
			F[0]=self.envelope_flux(T[0])
		# calculate dTdt
		dTdt=numpy.array([0.0]*self.ngrid)
		dTdt[:-1] = (self.eps[:-1] + (F[1:]-F[:-1])/(self.dx*self.P[:-1]))*self.g/cv[:-1]
		return dTdt
		
	def envelope_flux(self,T):
		return 10.0**self.TbTeff(log10(T))

	def crust_heating(self,P):
		# simple "smeared out" heating function, 1.2MeV in inner crust, 0.2MeV in outer crust
		heat = self.heat_source(P,1e16,1e17,1.7,self.dx)
		heat += self.heat_source(P,3e12,3e15,0.2,self.dx)
		# shallow heat source
		if self.shallow_heating:
			heat += self.heat_source(P,self.shallow_y/3.0,self.shallow_y*3.0,self.shallow_Q,self.dx)
		return heat

	def heat_source(self,double P,double y1,double y2,double Q, double dx):
		cdef double heat = 0.0
		cdef double P1 = P*exp(-0.5*dx)
		cdef double P2 = P*exp(0.5*dx)
		cdef double geff = 2.28e14  # gravity used to convert column to pressure
		if P1 > y1*geff and P2 < y2*geff:  # we are within the heating zone
			heat=8.8e4*Q*9.64e17/(P*log(y2/y1))
		if P1 < y1*geff and P2 < y2*geff and y1*geff<P2:   # left hand edge of heated region
			heat=8.8e4*Q*9.64e17/(P*log(y2/y1))
			heat *= log(P2/(y1*geff))/dx
		if P1 > y1*geff and P2 > y2*geff and y2*geff > P1: # right hand edge of heated region
			heat=8.8e4*Q*9.64e17/(P*log(y2/y1))
			heat *= log(y2*geff/P1)/dx
		return heat
		
	def	thickness(self):
		thickness = 0.0
		for item in self.grid:
			thickness+=item['P']*self.dx/item['rho']
		return thickness/self.g
