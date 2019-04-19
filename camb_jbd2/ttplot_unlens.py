from numpy import *
from scipy import *
from pylab import *
from astropy.io import fits
import matplotlib.pyplot as pyplot

grid1 = genfromtxt('test_scalCls.dat')
grid_lamda = genfromtxt('test_scalCls_lamda.dat')
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
tt_lamda = grid_lamda[:,1]

f,axarr = pyplot.subplots(2, sharex='True')
axarr[0].plot(l,tt_lamda,'k-',label = r"$\Lambda CDM$")
axarr[0].plot(l,tt_model_new,'r-',label = r"$\omega_{\rm{BD}}$")
axarr[1].plot(l,(tt_model_new-tt_lamda)/tt_lamda,'r-',label = r"$\omega_{\rm{BD}}$")
#axarr[0].plot(l_data_low, tt_data_low,'.b')
#axarr[0].plot(l_data_high, tt_data_high,'.b')
axarr[0].errorbar(l_data_low, tt_data_low,yerr=[err_data_up,err_data_low],fmt='k.')
axarr[0].errorbar(l_data_high, tt_data_high,err_data_high,fmt='k.')
#axarr[1].errorbar(l_data_low, tt_data_low-tt_data_low,yerr=[err_data_up,err_data_low],fmt='k.')
#axarr[1].errorbar(l_data_high, tt_data_high-tt_data_high,err_data_high,fmt='k.')
xlabel('l',fontsize=15)
axarr[0].set_ylabel(r'$\rm{l}(\rm{l}+1)C\rm{l}^{\rm{TT}}/2\pi$',fontsize=15)
axarr[1].set_ylabel(r'$\Delta_{\rm{rel}}$',fontsize=15)
#axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim([2.0,l[-1]])
axarr[0].set_xscale('log')
axarr[1].set_xscale('log')
axarr[0].set_ylim([0.0,6000.0])
axarr[0].legend(loc = 0)
axarr[1].legend(loc = 0)
savefig('cmb_ttcomp.pdf')
show()



