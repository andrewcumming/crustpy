from crust import *
import numpy
import matplotlib.pyplot as plt
import time

# These are the 1659 temperature measurements
t0 = 52159.5
tobs = numpy.array([52197.8,52563.2,52712.2,52768.9,53560.0,53576.7,54583.8,56113])-t0
Teffobs = numpy.array([121,85,77,73,58,54,56,48.8])

# Initialize the crust
crust = Crust(mass=1.6,radius=11.2,ngrid=30,Qimp=3.0,Tc=3e7)
print crust

# Set the top temperature and accrete
crust.set_top_temperature(4e8)
t, Teff = crust.evolve(time=2.5*365.0,mdot=0.1)
rho, TT = crust.temperature_profile()
K = numpy.array([item.Kcond() for item in crust.grid])
cve = numpy.array([item.CV_electrons() for item in crust.grid])
cvi = numpy.array([item.CV_ions() for item in crust.grid])
cv = numpy.array([item.CV() for item in crust.grid])

# Now turn off accretion and cool
t, Teff = crust.evolve(time=10000,mdot=0.0)
rho2, TT2 = crust.temperature_profile()

if 1:
	# Plot the temperature profile and lightcurve
	plt.subplot(2,2,1)
	plt.loglog(rho,TT)
	plt.loglog(rho2,TT2)

	ax=plt.subplot(2,2,2)
	plt.plot(t,Teff)
	plt.plot(tobs,Teffobs,'ro')
	plt.xlim([1.0,1e4])
	plt.ylim([50.0,130.0])
	ax.set_xscale('log')

	plt.subplot(2,2,3)
	plt.loglog(rho,K)

	plt.subplot(2,2,4)
	plt.loglog(rho,cvi,'k-')
	plt.loglog(rho,cve,'k--')
	plt.loglog(rho,cve+cvi,'k')
	plt.loglog(rho,cv,'b')

	plt.savefig('crustcool.pdf')
	print 'Saved plots to crustcool.pdf'



if 0:
	for item in crust.grid:
		print item, 'K=%g' % (item.Kcond()), 'CV=%g' % (item.CV())

	rho = numpy.array([item.rho for item in crust.grid])
	cve = numpy.array([item.CV_electrons() for item in crust.grid])
	cvi = numpy.array([item.CV_ions() for item in crust.grid])
	K = numpy.array([item.Kcond() for item in crust.grid])
	Kp = numpy.array([item.Kcond()*(item.fep()+item.feQ())/item.fep() for item in crust.grid])

	plt.subplot(2,1,1)
	plt.loglog(rho,K)
	plt.loglog(rho,Kp,'--')
	plt.subplot(2,1,2)
	plt.loglog(rho,cvi,'r')
	plt.loglog(rho,cve,'b')
	plt.loglog(rho,cve+cvi,'k')
	plt.show()
	