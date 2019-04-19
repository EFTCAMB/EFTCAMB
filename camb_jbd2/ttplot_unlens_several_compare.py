from numpy import *
from scipy import *
from pylab import *
from astropy.io import fits
import matplotlib.pyplot as pyplot

grid1 = genfromtxt('test_scalCls_75.dat')
grid2 = genfromtxt('test_scalCls_skordis_75.dat')
data = fits.open('COM_PowerSpect_CMB_R2.02.fits')
print data.info()
cmb_tt_low = data[1].data
cmb_tt_high = data[7].data

#Observational data
num = 27
num2 = 82
l_data_low = zeros(num,dtype=int)
l_data_high = zeros(num2,dtype=int)
tt_data_low = zeros(num,dtype = double)
tt_data_high = zeros(num2,dtype = double)
err_data_low = zeros(num,dtype = double)
err_data_up = zeros(num,dtype = double)
err_data_high = zeros(num2,dtype = double)

for i in range(0,27):
	l_data_low[i] = cmb_tt_low[i][0]
	tt_data_low[i] = cmb_tt_low[i][1]
	err_data_low[i] = array(cmb_tt_low[i][2]).T
	err_data_up[i] = array(cmb_tt_low[i][3]).T
for i in range(0,82):
	l_data_high[i] = cmb_tt_high[i][0]
	tt_data_high[i] = cmb_tt_high[i][3]
	err_data_high[i] = array(cmb_tt_high[i][4]).T
#theoretical data
l = grid1[:,0]
tt_model_new = grid1[:,1]
tt_model_new_2 = grid2[:,1]

figure(1)
plot(l,(tt_model_new),'r-',label = r"$\rm{EFTCAMB}$")
plot(l,(tt_model_new_2),'b-',label = r'$\rm{SKORDIS}$')
xlabel(r'$\ell$',fontsize=15)
ylabel(r'$C_{\ell}^{\rm{TT}}$',fontsize=15)
#axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim([2.0,l[-1]])
xscale('log')
legend(loc = 0)
savefig('cmb_ttcomp_1.pdf')

figure(2)
plot(l,(tt_model_new-tt_model_new_2)/tt_model_new_2*100,'r-',label = r"$\omega_{\rm{BD}} = 75$")
xlabel(r'$\ell$',fontsize=15)
ylabel(r'$\Delta C_{\ell}^{\rm{TT}}/C_{\ell}^{Avilez TT} (\%)$',fontsize=15)
#axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim([2.0,l[-1]])
xscale('log')
legend(loc = 0)
savefig('cmb_ttcomp_2.pdf')
show()



