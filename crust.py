import ns
import eos
import numpy

class Crust:
	def __init__(self,mass=1.4,radius=12.0,ngrid=30,Tc=1e8,Qimp=1):
		self.mass = mass
		self.radius = radius
		self.g = ns.grav(mass,radius)
		self.ngrid=ngrid
		self.Tc = Tc
		self.Qimp = Qimp
		P1 = 1e12*2.28e14
		P2 = 6.5e32
		self.grid = []
		for i in range(0,ngrid):
			P = numpy.log10(P1) + numpy.log10(P2/P1)*i/(ngrid-1)
			(A,Z,Yn) = self.composition(P)
			grid_point = eos.Eos(P=10.0**P,Yn=Yn,A=A,Z=Z,Qimp=Qimp)
			self.grid.append(grid_point)
						
	def __str__(self):
		description="M=%g R=%g g=%g ngrid=%d Tc=%g Qimp=%g" % (self.mass,self.radius,self.g,self.ngrid,self.Tc,self.Qimp)
		return description

	def composition(self,P):
		# composition as a function of pressure from HZ 1990
		Acella=[56.0,56.0,56.0,56.0,56.0,56.0,56.0,56.0,112.0,112.0,112.0,112.0,112.0,224.0,224.0,224.0,224.0,448.0,448.0]
		Aa=[56.0,56.0,56.0,56.0,56.0,52.0,46.0,40.0,68.0,62.0,56.0,50.0,44.0,66.0,60.0,54.0,48.0,96.0,88.0];
		Za=[26.0,24.0,22.0,20.0,18.0,16.0,14.0,12.0,20.0,18.0,16.0,14.0,12.0,18.0,16.0,14.0,12.0,24.0,22.0];
		# Maximum pressure at which that composition is present
		Pmaxa = [7.235e26,9.569e27,1.152e29,4.747e29,1.361e30,1.980e30,2.253e30,2.637e30,2.771e30,3.216e30,3.825e30,4.699e30,6.043e30,7.233e30,9.238e30,1.228e31,1.602e31,1.613e31,1e33]
		# Find the composition in the first entry in the table where P>Pmax.
		(Acell,A,Z) = [(Acell,A,Z) for (Acell,A,Z,Pmax) in zip(Acella,Aa,Za,Pmaxa) if Pmax>10.0**P][0]
		Yn = (Acell-A)/Acell
		return (A,Z,Yn)

	def evolve(self,time=10000.0,mdot=0.1):
		print "Evolving in time for %g days at mdot=%g" % (time,mdot)


		
if __name__ == '__main__':
	crust = Crust(1.4,12.0,ngrid=30,Qimp=1.0)
	print crust
	for item in crust.grid:
		print item, 'K=%g' % (item.Kcond()), 'CV=%g' % (item.CV())
