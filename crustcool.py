from crust import *
import numpy
import matplotlib.pyplot as plt

crust = Crust(mass=1.4,radius=12.0,ngrid=30,Qimp=1.0)

print crust
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