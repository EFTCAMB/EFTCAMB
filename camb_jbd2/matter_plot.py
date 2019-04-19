from numpy import *
from scipy import *
from pylab import *
from astropy.io import fits
import matplotlib.pyplot as pyplot


grid_matter_model = genfromtxt('test_matterpower_75.dat')
grid_matter_skordis = genfromtxt('test_matterpower_skordis_75.dat')
k = grid_matter_model[:,0]
k_2 = grid_matter_skordis[:,0]
#l = k[:]*9872.05*0.6781 #5979.31
matter_model = grid_matter_model[:,1]
matter_skordis = grid_matter_skordis[:,1]

figure(1)
plot(k,matter_model,label=r"$\rm{EFTCAMB}$")
plot(k_2,matter_skordis,label=r"$\rm{SKORDIS}$")
pyplot.xscale('log')
pyplot.yscale('log')
xlabel(r'k (h/Mpc)',fontsize=15)
xlim([1.0e-4,0.3])
ylim([7.0e2,3.0e4])
ylabel(r'$P(k)$',fontsize=15)
#ticklabel_format(style='sci', axis='y', scilimits=(0,0))
savefig('matter_comp.pdf')
legend(loc=0)
show()
