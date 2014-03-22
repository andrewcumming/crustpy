import ns
import eos
import numpy
from scipy.integrate import odeint
from scipy import interpolate

class Crust:
	def __init__(self,mass=1.4,radius=12.0,ngrid=30,Tc=1e8,Qimp=1):
		self.mass = mass
		self.radius = radius
		self.g,self.zz = ns.grav(mass,radius)
		self.ngrid=ngrid
		self.Tc = Tc
		self.Qimp = Qimp
		P1 = 1e12*2.28e14
		P2 = 6.5e32
		self.grid = []
		self.dx = numpy.log(P2/P1)/(ngrid-1)
		for i in range(0,ngrid):
			P = numpy.exp(numpy.log(P1) + self.dx*i)
			(A,Z,Yn) = self.composition(P)
			grid_point = eos.Eos(P=P,Yn=Yn,A=A,Z=Z,Qimp=Qimp,T=Tc)
			self.grid.append(grid_point)
		# Read in the envelope Tb-Teff relation
		F = numpy.array([float(line.split()[3]) for line in open('TbTeff','r')])
		T = numpy.array([float(line.split()[2]) for line in open('TbTeff','r')])
		self.TbTeff = interpolate.interp1d(T,F)
								
	def __str__(self):
		return "M=%g R=%g g=%g ngrid=%d Tc=%g Qimp=%g" % (self.mass,self.radius,self.g,self.ngrid,self.Tc,self.Qimp)

	def temperature_profile(self):
		return numpy.array([item.rho for item in self.grid]),numpy.array([item.T for item in self.grid])

	def set_top_temperature(self,T):
		self.grid[0].update_T(T)

	def set_core_temperature(self,T):
		self.grid[-1].update_T(T)

	def composition(self,P):
		# composition as a function of pressure from HZ 1990
		Acella=[56.0,56.0,56.0,56.0,56.0,56.0,56.0,56.0,112.0,112.0,112.0,112.0,112.0,224.0,224.0,224.0,224.0,448.0,448.0]
		Aa=[56.0,56.0,56.0,56.0,56.0,52.0,46.0,40.0,68.0,62.0,56.0,50.0,44.0,66.0,60.0,54.0,48.0,96.0,88.0];
		Za=[26.0,24.0,22.0,20.0,18.0,16.0,14.0,12.0,20.0,18.0,16.0,14.0,12.0,18.0,16.0,14.0,12.0,24.0,22.0];
		# Maximum pressure at which that composition is present
		Pmaxa = [7.235e26,9.569e27,1.152e29,4.747e29,1.361e30,1.980e30,2.253e30,2.637e30,2.771e30,3.216e30,3.825e30,4.699e30,6.043e30,7.233e30,9.238e30,1.228e31,1.602e31,1.613e31,1e33]
		# Find the composition in the first entry in the table where P>Pmax.
		(Acell,A,Z) = [(Acell,A,Z) for (Acell,A,Z,Pmax) in zip(Acella,Aa,Za,Pmaxa) if Pmax>P][0]
		Yn = (Acell-A)/Acell
		return (A,Z,Yn)

	def evolve(self,time=10000.0,mdot=0.1,n=100):
		print "Evolving in time for %g days at mdot=%g" % (time,mdot)
		inic = numpy.array([item.T for item in self.grid])
		times = 10.0**(numpy.arange(n)*numpy.log10(time)/(n-1))
		result = odeint(self.derivs,inic,times*3600.0*24.0,rtol=1e-6,atol=1e-6,args=(mdot,))
		for i,T in enumerate(result[-1,:]):
			self.grid[i].update_T(T)
		Teff = [(self.envelope_flux(T)/5.67e-5)**0.25*1.38e-16/(1.6e-12*self.zz) for T in result[:,0]]
		# returns the redshifted lightcurve
		return times*self.zz, Teff
		
	def derivs(self,T,time,mdot):
		# calculate fluxes
		g=self.g*1e14
		K = [item.rho*g*item.Kcond()/item.P for item in self.grid]
		# F[j] is the flux at j+1/2
		F = [ 0.5*(K[j]+K[j+1])*(T[j+1]-T[j])/self.dx for j in range(0,self.ngrid-1)]
		eps = [mdot * g* self.crust_heating(item.P) for item in self.grid]
		# calculate dTdt
		dTdt = [(eps[j] + g*(F[j]-F[j-1])/(self.dx*self.grid[j].P))/self.grid[j].CV() for j in range(1,self.ngrid-1)]
		if mdot>0.0:
			dTdt0=0.0
		else:
			dTdt0 = g*(F[0]-self.envelope_flux(T[0]))/(self.dx*self.grid[0].P*self.grid[0].CV())
		return [dTdt0]+dTdt+[0.0]

	def envelope_flux(self,T):
		return 10.0**self.TbTeff(numpy.log10(T))

	def crust_heating(self,P):
		# simple "smeared out" heating function, 1.2MeV in inner crust, 0.2MeV in outer crust
		if P>=1e16*2.28e14 and P<=1e17*2.28e14:
			return 8.8e4*1.7*9.64e17/(P*numpy.log(1e17/1e16));
	 	if P>=3e12*2.28e14 and P<3e15*2.28e14:
			return 8.8e4*0.2*9.64e17/(P*numpy.log(3e15/3e12));
		return 0.0


if __name__ == '__main__':
	crust = Crust(1.4,12.0,ngrid=30,Qimp=1.0)
	print crust
	for item in crust.grid:
		print item, 'K=%g' % (item.Kcond()), 'CV=%g' % (item.CV())
	crust.derivs([1.0],1.0)