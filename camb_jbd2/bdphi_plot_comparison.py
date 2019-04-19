from scipy import *
from numpy import *
from pylab import *
import matplotlib.pyplot as pyplot

a = genfromtxt('a_evol.txt')
phi = genfromtxt('phi_a.txt')
hubble = genfromtxt('hubble.txt')
a_2 = genfromtxt('a_evo.txt')
phi_2 = genfromtxt('phi_evo.txt')
hubble_2 = genfromtxt('hubble_evo.txt')

figure(1)
plot(a,phi,label = r'$\rm{EFTCAMB}$')
plot(a_2,phi_2,label = r'$\rm{SKORDIS}$')
pyplot.xscale('log')
xlim([1.0e-10,1.0])
xlabel('a',fontsize=15)
ylabel(r'$\phi$',fontsize=15)
legend(loc=0)
savefig('phi_comp.pdf')

figure(2)
plot(a,hubble,label = r'$\rm{EFTCAMB}$')#-hubble_lamda)/hubble_lamda)
plot(a_2,hubble_2,label = r'$\rm{SKORDIS}$')
pyplot.xscale('log')
#pyplot.yscale('log')
xlim([1.0e-1,1.0])
ylim([60,150.0])
xlabel('a',fontsize=15)
ylabel(r'$H(a)$',fontsize=15)
legend(loc=0)
savefig('Hubble_comp.pdf')

show()
