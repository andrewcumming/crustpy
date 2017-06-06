from crust import *
import numpy
import matplotlib.pyplot as plt
import matplotlib
import time
matplotlib.rcParams.update({'font.size': 8})

# These are the 1659 temperature measurements
t0 = 52159.5
tobs = numpy.array([52197.8,52563.2,52712.2,52768.9,53560.0,53576.7,54583.8,56113])-t0
Teffobs = numpy.array([121,85,77,73,58,54,56,48.8])

# Initialize the crust
crust = Crust(mass=1.62,radius=11.2,ngrid=50,Qimp=6.0,Tc=3.1e7)
print(crust)

# Set the top temperature and accrete
crust.set_top_temperature(4.7e8)
t, Teff = crust.evolve(time=2.5*365.0,mdot=0.1)

# store the initial temperature, conductivity and CV profile for the plots
rho, TT = crust.temperature_profile()
K = numpy.array([item['Kcond'] for item in crust.grid])
cve = numpy.array([item['CV_electrons'] for item in crust.grid])
cvi = numpy.array([item['CV_ions'] for item in crust.grid])
cv = numpy.array([item['CV'] for item in crust.grid])
# print out crust profile just before cooling begins:
#for item in crust.grid:
#	print(item, 'K=%g' % (item['Kcond']), 'CV=%g' % (item['CV']))

# Now turn off accretion and cool
t, Teff = crust.evolve(time=10000,mdot=0.0)

# Final temperature profile to plot
rho2, TT2 = crust.temperature_profile()

if 1:
	# Plot the temperature profile and lightcurve
	plt.subplot(2,2,1)
	plt.xlim([1e8,2e14])
	plt.loglog(rho,TT)
	plt.loglog(rho2,TT2)
	plt.xlabel(r'$\rho (\mathrm{g/cm^3})$')
	plt.ylabel(r'$T (\mathrm{K})$')

	ax=plt.subplot(2,2,2)
	plt.plot(t,Teff)
	plt.plot(tobs,Teffobs,'ro',markersize=2)
	plt.xlim([10.0,1e4])
	plt.ylim([50.0,130.0])
	ax.set_xscale('log')
	plt.xlabel(r'$t (\mathrm{d})$')
	plt.ylabel(r'$T_{eff} (\mathrm{eV})$')

	plt.subplot(2,2,3)
	plt.xlim([1e8,2e14])
	plt.loglog(rho,K)
	plt.xlabel(r'$\rho (\mathrm{g/cm^3})$')
	plt.ylabel(r'$K (\mathrm{cgs})$')

	plt.subplot(2,2,4)
	plt.xlim([1e8,2e14])
	plt.loglog(rho,cv)
	plt.loglog(rho,cvi,':')
	plt.loglog(rho,cve,'--')
	#plt.loglog(rho,cve+cvi)
	plt.xlabel(r'$\rho (\mathrm{g/cm^3})$')
	plt.ylabel(r'$C_V (\mathrm{erg/g/K})$')

	plt.tight_layout()

	#plt.show()
	plt.savefig('crustcool.pdf')
	print('Saved plots to crustcool.pdf')
