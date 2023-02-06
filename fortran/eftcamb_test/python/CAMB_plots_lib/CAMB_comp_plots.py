# ***************************************************************************************

from .CMB_plots        import CMB_plots
from scipy.interpolate import UnivariateSpline

import numpy             as np
import matplotlib.pyplot as plt
import math
import os

# ***************************************************************************************

""" Hard coded options """

cutoff = 10.0**(-16)

# ***************************************************************************************

class CAMB_results_compare_plot:
    """
    class that contains the necessary tools to plot the comparison of two CAMB runs
    """

    # general plot settings:
    color1      = 'red'                # default color of the line of the first model
    color2      = 'blue'               # default color of the line of the second model
    color_compa = 'green'              # default color of the line of the difference
    x_size_inc  = 8.30                 # x size of the final plots, in inches
    y_size_inc  = 11.7                 # y size of the final plots, in inches

    def __init__(self, root1, root2, outpath,
                 name1='', name2=''):
        """
        Class constructor:
            root1     = name of the first CAMB run
            root2     = name of the second CAMB run
            outpath   = path to the output folder
            name1     = (optional) specifies the name of the first model. Used for the legend.
            name2     = (optional) specifies the name of the second model. Used for the legend.
        """

        # store the constructor options:
        self.root1    = root1
        self.root2    = root2
        self.outpath  = outpath

        # extract the name of the models from the roots:
        self.name1 = root1.split("/")[len(root1.split("/"))-1]
        self.name2 = root2.split("/")[len(root2.split("/"))-1]

        # set the human readable name
        if name1=='':
            self.name_h1 = ''.join(i for i in self.name1.replace('_',' ') if not i.isdigit())
        else:
            self.name_h1 = name1

        if name2=='':
            self.name_h2 = ''.join(i for i in self.name2.replace('_',' ') if not i.isdigit())
        else:
            self.name_h2 = name2

        # load the data:
        self.scalCls_1           = np.loadtxt(root1+'_scalCls.dat')
        self.scalarCovCls_1        = np.loadtxt(root1+'_scalarCovCls.dat')

        self.scalCls_2           = np.loadtxt(root2+'_scalCls.dat')
        self.scalarCovCls_2        = np.loadtxt(root2+'_scalarCovCls.dat')
        
        # get the lensing spectra:
        if os.path.exists(root1+'_lensedCls.dat') and os.path.exists(root2+'_lensedCls.dat'):
            self.lensing = True
            self.lensedCls_1         = np.loadtxt(root1+'_lensedCls.dat')
            self.lenspotentialCls_1  = np.loadtxt(root1+'_lenspotentialCls.dat')
            self.lensedCls_2         = np.loadtxt(root2+'_lensedCls.dat')
            self.lenspotentialCls_2  = np.loadtxt(root2+'_lenspotentialCls.dat')
        else:
            self.lensing = False

        if os.path.exists(root1+'_matterpower.dat') and os.path.exists(root2+'_matterpower.dat'):
            self.transfer = True
            self.matterpower_1       = np.loadtxt(root1+'_matterpower.dat')
            self.transfer_func_1     = np.loadtxt(root1+'_transfer_out.dat')
            self.matterpower_2       = np.loadtxt(root2+'_matterpower.dat')
            self.transfer_func_2     = np.loadtxt(root2+'_transfer_out.dat')
            # try to get the labels of the transfer function:
            self.transfer_labels = None
            infile = open(root1+'_transfer_out.dat', 'r')
            firstLine = infile.readline()
            if firstLine[0]=='#':
                firstLine = firstLine.split(' ')
                firstLine = list(filter(None,list(firstLine)))
                firstLine.pop(0)
                firstLine.pop(0)
                firstLine.pop(-1)
                self.transfer_labels = firstLine
        else:
            self.transfer = False

        if os.path.exists(root1+'_tensCls.dat') and os.path.exists(root2+'_tensCls.dat'):
            self.tensor = True
            self.tensCls_1           = np.loadtxt(root1+'_tensCls.dat')
            self.totCls_1            = np.loadtxt(root1+'_totCls.dat')
            self.tensCls_2           = np.loadtxt(root2+'_tensCls.dat')
            self.totCls_2            = np.loadtxt(root2+'_totCls.dat')
        else:
            self.tensor = False

        if self.lensing and self.tensor:
            self.lensedtotCls_1      = np.loadtxt(root1+'_lensedtotCls.dat')
            self.lensedtotCls_2      = np.loadtxt(root2+'_lensedtotCls.dat')

    def plot_compare_scalCls(self):
        """
        Plots and saves the comparison of all the scalar Cls in a unique image
        """

        # number of Cls in the two files:
        num1     = self.scalCls_1.shape[1]-1
        num2     = self.scalCls_2.shape[1]-1
        # protection against different runs
        if num1!=num2:
            print( 'wrong number of Cls' )
            return
        # add the matter power spectrum if required:
        if self.transfer: num1 += 1
        # values of l:
        xvalues_1 = self.scalCls_1[:,0]
        xvalues_2 = self.scalCls_2[:,0]
        # get minimuma and maximum l values:
        l_min = int(np.amax( [np.amin(xvalues_1), np.amin(xvalues_2)] ))
        l_max = int(np.amin( [np.amax(xvalues_1), np.amax(xvalues_2)] ))
        # get the l values:
        l_values = np.linspace(l_min, l_max, l_max-l_min)

        # set up the plots:
        plots_1       = CMB_plots()
        plots_2       = CMB_plots()
        plots_compa   = CMB_plots()

        plots_1.color          = self.color1
        plots_2.color          = self.color2
        plots_compa.color      = self.color_compa
        plots_compa.comparison = True
        plots_compa.axes_label_position = 'right'
        fig = plt.gcf()

        # do the plots:
        for ind in range(1,num1+1):

            # distribute the plots in the figure:
            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            # do the l sampling and comparison:
            if not ( ind == num1 and self.transfer ):
                # get the y values:
                yvalues_1    = self.scalCls_1[:,ind]
                yvalues_2    = self.scalCls_2[:,ind]
                # interpolate the two spectra:
                spline_1 = UnivariateSpline( xvalues_1, yvalues_1, s=0)
                spline_2 = UnivariateSpline( xvalues_2, yvalues_2, s=0)
                # evaluate on the l grid:
                yvalues_1 = spline_1(l_values)
                yvalues_2 = spline_2(l_values)
                # protect against zeroes:
                yvalues_1_temp = np.abs( yvalues_1 )
                yvalues_2_temp = np.abs( yvalues_2 )
                try:
                    min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                    min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
                except:
                    min2val_1      = cutoff
                    min2val_2      = cutoff
                np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
                np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
                # computation of comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)
                # account for variance of TE:
                if ind == 3:
                    yvalues_comp = (yvalues_1 - yvalues_2)/np.sqrt(2./(2.*l_values+1.))
                    _temp_ClTT = UnivariateSpline( xvalues_1, self.scalCls_1[:,1], s=0)(l_values)
                    _temp_ClEE = UnivariateSpline( xvalues_1, self.scalCls_1[:,2], s=0)(l_values)
                    yvalues_comp = yvalues_comp/np.sqrt( _temp_ClTT*_temp_ClEE+yvalues_1**2.)
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<cutoff, [cutoff])
            else: # matter power spectrum case
                # get the k values:
                xvalues_1      = self.matterpower_1[:,0]
                xvalues_2      = self.matterpower_2[:,0]
                # get the y values:
                yvalues_1    = self.matterpower_1[:,1]
                yvalues_2    = self.matterpower_2[:,1]
                # get the k grid:
                k_min = np.amax( [np.amin(xvalues_1), np.amin(xvalues_2)] )
                k_max = np.amin( [np.amax(xvalues_1), np.amax(xvalues_2)] )
                k_values = np.logspace(np.log10(k_min), np.log10(k_max), 1000)
                # interpolate the two spectra:
                spline_1 = UnivariateSpline( xvalues_1, yvalues_1, s=0)
                spline_2 = UnivariateSpline( xvalues_2, yvalues_2, s=0)
                # evaluate on the l grid:
                yvalues_1 = spline_1(k_values)
                yvalues_2 = spline_2(k_values)
                # protect against zeroes:
                yvalues_1_temp = np.abs( yvalues_1 )
                yvalues_2_temp = np.abs( yvalues_2 )
                try:
                    min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                    min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
                except:
                    min2val_1      = cutoff
                    min2val_2      = cutoff
                np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
                np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<cutoff, [cutoff])

            # make the plots:
            if ind == 1: # TT power spectrum:

                plots_1.TT_plot(temp, l_values, yvalues_1)
                plots_2.TT_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.TT_plot(temp_comp, l_values, yvalues_comp)
                temp_comp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2: # EE power spectrum:

                plots_1.EE_plot(temp, l_values, yvalues_1)
                plots_2.EE_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.EE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3: # TE power spectrum:

                plots_1.TE_plot(temp, l_values, yvalues_1)
                plots_2.TE_plot(temp, l_values, yvalues_2)
                plots_compa.TE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('TE power spectrum')

            elif ind == 4 and self.lensing: # CMB lensing power spectrum:

                plots_1.Phi_plot(temp, l_values, yvalues_1)
                plots_2.Phi_plot(temp, l_values, yvalues_2)
                plots_compa.Phi_plot(temp_comp, l_values, yvalues_comp*100.)
                temp.set_title('$\phi$ power spectrum')

            elif ind == 5 and self.lensing: # CMB lensing - Temperature power spectrum:

                plots_1.PhiT_plot(temp, l_values, yvalues_1)
                plots_2.PhiT_plot(temp, l_values, yvalues_2)
                plots_compa.PhiT_plot(temp_comp, l_values, yvalues_comp*100.)
                temp.set_title('$\phi$T power spectrum')

            elif ind == num1 and self.transfer: # matter power spectrum:

                plots_1.Matter_plot(temp, k_values, yvalues_1)
                plots_2.Matter_plot(temp, k_values, yvalues_2)
                plots_compa.Matter_plot(temp_comp, k_values, yvalues_comp)
                temp.set_title('Matter power spectrum')

            else: # generic Cl comparison:

                plots_1.Generic_Cl(temp, l_values, yvalues_1)
                plots_2.Generic_Cl(temp, l_values, yvalues_2)
                plots_compa.Generic_Cl(temp_comp, l_values, yvalues_comp)

        # set the size of the image
        fig.set_size_inches( self.x_size_inc, self.y_size_inc)
        # set a tight layout
        fig.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)
        # set the global title
        plt.suptitle(self.name_h1+' VS '+self.name_h2+
                     ' comparison of scalar Cls', fontsize=16)
        # set the global legend
        fig.legend( handles = [plots_1.TT_p, plots_2.TT_p, plots_compa.CV_p],
                    labels  = [self.name_h1, self.name_h2, 'Cosmic variance'],
                    loc='lower center', ncol=3 ,fancybox=True, edgecolor='k')
        # adjust the size of the plot
        fig.subplots_adjust(top=0.92, bottom=0.08)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_scalCls_comp.pdf')
        plt.clf()
        plt.close("all")

    def plot_compare_scalarCovCls(self):
        """
        Plots and saves the comparison of all the scalar Cov Cls in a unique image
        """

        # number of Cls:
        num1     = self.scalarCovCls_1.shape[1]-1
        num2     = self.scalarCovCls_2.shape[1]-1
        # protection against different runs
        if num1!=num2:
            print( 'wrong number of Cls' )
            return

        # size of the Cl Cov matrix:
        num1     = int(math.sqrt(num1))
        num_tot  = sum(range(1,num1+1))

        # set up the plots:
        plots_1       = CMB_plots()
        plots_2       = CMB_plots()
        plots_compa   = CMB_plots()

        plots_1.color          = self.color1
        plots_2.color          = self.color2
        plots_compa.color      = self.color_compa
        plots_compa.comparison = True
        plots_compa.axes_label_position = 'right'
        fig = plt.gcf()

        # setup a dictionary with the names of the Cls
        dict = { 1: 'T', 2: 'E', 3: '$\phi$'}
        for i in range(4, num1+1):
            dict[i] = 'W'+str(i)

        # values of the multipoles:
        xvalues_1 = self.scalarCovCls_1[:,0]
        xvalues_2 = self.scalarCovCls_2[:,0]
        # get the l values:
        l_min = int( np.amax( [np.amin(xvalues_1), np.amin(xvalues_2)] ) )
        l_max = int( np.amin( [np.amax(xvalues_1), np.amax(xvalues_2)] ) )
        # get the l values:
        l_values = np.linspace(l_min, l_max, l_max-l_min)

        # other stuff:
        ind_tot = 0
        # do the plots:
        for ind in range(1,num1+1):
            for ind2 in range(1, ind+1):

                ind_tot += 1
                # place the plots:
                temp      = plt.subplot2grid((num_tot,2), (ind_tot-1, 0))
                temp_comp = plt.subplot2grid((num_tot,2), (ind_tot-1, 1))
                # values of the Cls:
                col          = ind + num1*(ind2-1)
                # get the y values:
                yvalues_1    = self.scalarCovCls_1[:,col]
                yvalues_2    = self.scalarCovCls_2[:,col]
                # interpolate the two spectra:
                spline_1 = UnivariateSpline( xvalues_1, yvalues_1, s=0)
                spline_2 = UnivariateSpline( xvalues_2, yvalues_2, s=0)
                # evaluate on the l grid:
                yvalues_1 = spline_1(l_values)
                yvalues_2 = spline_2(l_values)
                # protect against zeroes:
                yvalues_1_temp = np.abs( yvalues_1 )
                yvalues_2_temp = np.abs( yvalues_2 )
                try:
                    min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                    min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
                except:
                    min2val_1      = cutoff
                    min2val_2      = cutoff
                np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
                np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<cutoff, [cutoff])

                # make the plots:
                plots_1.Generic_Cl(temp, l_values, yvalues_1)
                plots_2.Generic_Cl(temp, l_values, yvalues_2)
                plots_compa.Generic_Cl(temp_comp, l_values, yvalues_comp)

                temp.set_title(dict[ind2]+dict[ind]+' power spectrum')

        # set the size of the image
        fig.set_size_inches( self.x_size_inc, self.y_size_inc/5.0*num_tot)
        # set a tight layout
        fig.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)

        # set the global legend
        fig.legend( handles = [plots_1.Generic_Cl_plot, plots_2.Generic_Cl_plot, plots_compa.CV_p],
                    labels  = [self.name_h1, self.name_h2, 'Cosmic variance'],
                    loc='lower center', ncol=3 ,fancybox=True, edgecolor='k')

        # set the global title
        plt.suptitle(self.name_h1+' VS '+self.name_h2+
                     ' comparison of scalar Cov Cls', fontsize=16)

        # adjust the size of the plot
        fig.subplots_adjust(top=0.92, bottom=0.08)
        #        fig.subplots_adjust(top=0.96, bottom=0.01)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_scalarCovCls_comp.pdf')
        plt.clf()
        plt.close("all")

    def plot_compare_lensedCls(self):
        """
        Plots and saves the comparison of all the lensed Cls in a unique image
        """

        # protection from direct calls if lensing is not included
        if not self.lensing: return

        # number of Cls:
        num1     = self.lensedCls_1.shape[1]-1
        num2     = self.lensedCls_2.shape[1]-1
        # protection against different runs
        if num1!=num2:
            print( 'wrong number of Cls' )
            return

        # get the x values:
        xvalues_1 = self.lensedCls_1[:,0]
        xvalues_2 = self.lensedCls_2[:,0]
        # get the l values:
        l_min = int(np.amax( [np.amin(xvalues_1), np.amin(xvalues_2)] ))
        l_max = int(np.amin( [np.amax(xvalues_1), np.amax(xvalues_2)] ))
        # get the l values:
        l_values = np.linspace(l_min, l_max, l_max-l_min)

        # set up the plots:
        plots_1       = CMB_plots()
        plots_2       = CMB_plots()
        plots_compa   = CMB_plots()

        plots_1.color          = self.color1
        plots_2.color          = self.color2
        plots_compa.color      = self.color_compa
        plots_compa.comparison = True
        plots_compa.axes_label_position = 'right'
        fig = plt.gcf()

        # do the plots:
        for ind in range(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            # get the y values:
            yvalues_1    = self.lensedCls_1[:,ind]
            yvalues_2    = self.lensedCls_2[:,ind]
            # interpolate the two spectra:
            spline_1 = UnivariateSpline( xvalues_1, yvalues_1, s=0)
            spline_2 = UnivariateSpline( xvalues_2, yvalues_2, s=0)
            # evaluate on the l grid:
            yvalues_1 = spline_1(l_values)
            yvalues_2 = spline_2(l_values)
            # protect against zeroes:
            yvalues_1_temp = np.abs( yvalues_1 )
            yvalues_2_temp = np.abs( yvalues_2 )
            try:
                min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
            except:
                min2val_1      = cutoff
                min2val_2      = cutoff
            np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
            np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
            # computation of the percentual comparison:
            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)
            # account for variance of TE:
            if ind == 4:
                yvalues_comp = (yvalues_1 - yvalues_2)/np.sqrt(2./(2.*l_values+1.))
                _temp_ClTT = UnivariateSpline( xvalues_1, self.lensedCls_1[:,1], s=0)(l_values)
                _temp_ClEE = UnivariateSpline( xvalues_1, self.lensedCls_1[:,2], s=0)(l_values)
                yvalues_comp = yvalues_comp/np.sqrt( _temp_ClTT*_temp_ClEE+yvalues_1**2.)
            # protection against values too small:
            np.place(yvalues_comp, abs(yvalues_comp)<cutoff, [cutoff])

            if ind == 1:

                plots_1.TT_plot(temp, l_values, yvalues_1)
                plots_2.TT_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.TT_plot(temp_comp, l_values, yvalues_comp)
                temp_comp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2:

                plots_1.EE_plot(temp, l_values, yvalues_1)
                plots_2.EE_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.EE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3:

                plots_1.BB_plot(temp, l_values, yvalues_1)
                plots_2.BB_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.BB_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('BB power spectrum')

            elif ind == 4:

                plots_1.TE_plot(temp, l_values, yvalues_1)
                plots_2.TE_plot(temp, l_values, yvalues_2)
                plots_compa.TE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('TE power spectrum')

            else:

                plots_1.TT_plot(temp, l_values, yvalues_1)
                plots_2.TT_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.TT_plot(temp_comp, l_values, yvalues_comp)

        # set the size of the image
        fig.set_size_inches( self.x_size_inc, self.y_size_inc)
        # set a tight layout
        fig.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)
        # set the global title
        plt.suptitle(self.name_h1+' VS '+self.name_h2+
                     ' comparison of lensed Cls', fontsize=16)
        # set the global legend
        fig.legend( handles = [plots_1.TT_p, plots_2.TT_p, plots_compa.CV_p],
                    labels  = [self.name_h1, self.name_h2, 'Cosmic variance'],
                    loc='lower center', ncol=3 ,fancybox=True, edgecolor='k')
        # adjust the size of the plot
        fig.subplots_adjust(top=0.92, bottom=0.08)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_lensedCls_comp.pdf')
        plt.clf()
        plt.close("all")

    def plot_compare_tensCls(self):
        """
        Plots and saves the comparison of all the tensor Cls in a unique image
        """

        # protection from direct calls if tensors are not included
        if not self.tensor: return

        # number of Cls:
        num1     = self.tensCls_1.shape[1]-1
        num2     = self.tensCls_2.shape[1]-1
        # protection against different runs
        if num1!=num2:
            print( 'wrong number of Cls' )
            return

        # get the x values:
        xvalues_1 = self.tensCls_1[:,0]
        xvalues_2 = self.tensCls_2[:,0]
        # get the l values:
        l_min = int(np.amax( [np.amin(xvalues_1), np.amin(xvalues_2)] ))
        l_max = int(np.amin( [np.amax(xvalues_1), np.amax(xvalues_2)] ))
        # get the l values:
        l_values = np.linspace(l_min, l_max, l_max-l_min)

        # set up the plots:
        plots_1       = CMB_plots()
        plots_2       = CMB_plots()
        plots_compa   = CMB_plots()

        plots_1.color          = self.color1
        plots_2.color          = self.color2
        plots_compa.color      = self.color_compa
        plots_compa.comparison = True
        plots_compa.axes_label_position = 'right'
        fig = plt.gcf()

        # do the plots:
        for ind in range(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            # get the y values:
            yvalues_1    = self.tensCls_1[:,ind]
            yvalues_2    = self.tensCls_2[:,ind]
            # interpolate the two spectra:
            spline_1 = UnivariateSpline( xvalues_1, yvalues_1, s=0)
            spline_2 = UnivariateSpline( xvalues_2, yvalues_2, s=0)
            # evaluate on the l grid:
            yvalues_1 = spline_1(l_values)
            yvalues_2 = spline_2(l_values)
            # protect against zeroes:
            yvalues_1_temp = np.abs( yvalues_1 )
            yvalues_2_temp = np.abs( yvalues_2 )
            try:
                min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
            except:
                min2val_1      = cutoff
                min2val_2      = cutoff
            np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
            np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
            # computation of the percentual comparison:
            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)
            # account for variance of TE:
            if ind == 4:
                yvalues_comp = (yvalues_1 - yvalues_2)/np.sqrt(2./(2.*l_values+1.))
                _temp_ClTT = UnivariateSpline( xvalues_1, self.tensCls_1[:,1], s=0)(l_values)
                _temp_ClEE = UnivariateSpline( xvalues_1, self.tensCls_1[:,2], s=0)(l_values)
                yvalues_comp = yvalues_comp/np.sqrt( _temp_ClTT*_temp_ClEE+yvalues_1**2.)
            # protection against values too small:
            np.place(yvalues_comp, abs(yvalues_comp)<cutoff, [cutoff])

            if ind == 1:

                plots_1.TT_plot(temp, l_values, yvalues_1)
                plots_2.TT_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.TT_plot(temp_comp, l_values, yvalues_comp)
                temp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2:

                plots_1.EE_plot(temp, l_values, yvalues_1)
                plots_2.EE_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.EE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3:

                plots_1.BB_plot(temp, l_values, yvalues_1)
                plots_2.BB_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.BB_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('BB power spectrum')

            elif ind == 4:

                plots_1.TE_plot(temp, l_values, yvalues_1)
                plots_2.TE_plot(temp, l_values, yvalues_2)
                plots_compa.TE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('TE power spectrum')

            else:

                plots_1.TT_plot(temp, l_values, yvalues_1)
                plots_2.TT_plot(temp, l_values, yvalues_2)
                plots_compa.TT_plot(temp_comp, l_values, yvalues_comp)

        # set the size of the image
        fig.set_size_inches( self.x_size_inc, self.y_size_inc)
        # set a tight layout
        fig.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)
        # set the global title
        plt.suptitle(self.name_h1+' VS '+self.name_h2+
                     ' comparison of tensor Cls', fontsize=16)
        # set the global legend
        fig.legend( handles = [plots_1.TT_p, plots_2.TT_p, plots_compa.CV_p],
                    labels  = [self.name_h1, self.name_h2, 'Cosmic variance'],
                    loc='lower center', ncol=3 ,fancybox=True, edgecolor='k')
        # adjust the size of the plot
        fig.subplots_adjust(top=0.92, bottom=0.08)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_tensCls_comp.pdf')
        plt.clf()
        plt.close("all")

    def plot_compare_totalCls(self):
        """
        Plots and saves the comparison of all the total (scalar + tensor) Cls in a unique image
        If lensing is included lensed Cls are used.
        """

        # protection from direct calls if tensors are not included
        if not self.tensor: return

        # decide what data to use:
        if self.lensing:
            data1 = self.lensedtotCls_1
            data2 = self.lensedtotCls_2
        else:
            data1 = self.totCls_1
            data2 = self.totCls_2

        # number of Cls:
        num1     = data1.shape[1]-1
        num2     = data2.shape[1]-1
        # protection against different runs
        if num1!=num2:
            print( 'wrong number of Cls' )
            return

        # get the x values:
        xvalues_1 = data1[:,0]
        xvalues_2 = data2[:,0]
        # get the l values:
        l_min = int( np.amax( [np.amin(xvalues_1), np.amin(xvalues_2)] ) )
        l_max = int( np.amin( [np.amax(xvalues_1), np.amax(xvalues_2)] ) )
        # get the l values:
        l_values = np.linspace(l_min, l_max, l_max-l_min)

        # set up the plots:
        plots_1       = CMB_plots()
        plots_2       = CMB_plots()
        plots_compa   = CMB_plots()

        plots_1.color          = self.color1
        plots_2.color          = self.color2
        plots_compa.color      = self.color_compa
        plots_compa.comparison = True
        plots_compa.axes_label_position = 'right'
        fig = plt.gcf()

        # do the plots:
        for ind in range(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            # get the y values:
            yvalues_1    = data1[:,ind]
            yvalues_2    = data2[:,ind]
            # interpolate the two spectra:
            spline_1 = UnivariateSpline( xvalues_1, yvalues_1, s=0)
            spline_2 = UnivariateSpline( xvalues_2, yvalues_2, s=0)
            # evaluate on the l grid:
            yvalues_1 = spline_1(l_values)
            yvalues_2 = spline_2(l_values)
            # protect against zeroes:
            yvalues_1_temp = np.abs( yvalues_1 )
            yvalues_2_temp = np.abs( yvalues_2 )
            try:
                min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
            except:
                min2val_1      = cutoff
                min2val_2      = cutoff
            np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
            np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
            # computation of the percentual comparison:
            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)
            # account for variance of TE:
            if ind == 4:
                yvalues_comp = (yvalues_1 - yvalues_2)/np.sqrt(2./(2.*l_values+1.))
                _temp_ClTT = UnivariateSpline( xvalues_1, data1[:,1], s=0)(l_values)
                _temp_ClEE = UnivariateSpline( xvalues_1, data1[:,2], s=0)(l_values)
                yvalues_comp = yvalues_comp/np.sqrt( _temp_ClTT*_temp_ClEE+yvalues_1**2.)
            # protection against values too small:
            np.place(yvalues_comp, abs(yvalues_comp)<cutoff, [cutoff])

            if ind == 1:

                plots_1.TT_plot(temp, l_values, yvalues_1)
                plots_2.TT_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.TT_plot(temp_comp, l_values, yvalues_comp)
                temp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2:

                plots_1.EE_plot(temp, l_values, yvalues_1)
                plots_2.EE_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.EE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3:

                plots_1.BB_plot(temp, l_values, yvalues_1)
                plots_2.BB_plot(temp, l_values, yvalues_2)
                yvalues_comp = yvalues_comp/np.sqrt(2./(2.*l_values+1.))
                plots_compa.BB_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('BB power spectrum')

            elif ind == 4:

                plots_1.TE_plot(temp, l_values, yvalues_1)
                plots_2.TE_plot(temp, l_values, yvalues_2)
                plots_compa.TE_plot(temp_comp, l_values, yvalues_comp)
                temp.set_title('TE power spectrum')

            else:

                plots_1.TT_plot(temp, l_values, yvalues_1)
                plots_2.TT_plot(temp, l_values, yvalues_2)
                plots_compa.TT_plot(temp_comp, l_values, yvalues_comp)


        # set the size of the image
        fig.set_size_inches( self.x_size_inc, self.y_size_inc)
        # set a tight layout
        fig.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)
        # set the global title
        if self.lensing:
            plt.suptitle(self.name_h1+' VS '+self.name_h2+' comparison of total lensed Cls', fontsize=16)
        else:
            plt.suptitle(self.name_h1+' VS '+self.name_h2+' comparison of total Cls', fontsize=16)
        # set the global legend
        fig.legend( handles = [plots_1.TT_p, plots_2.TT_p, plots_compa.CV_p],
                    labels  = [self.name_h1, self.name_h2, 'Cosmic variance'],
                    loc='lower center', ncol=3 ,fancybox=True, edgecolor='k')
        # adjust the size of the plot
        fig.subplots_adjust(top=0.92, bottom=0.08)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_totCls_comp.pdf')
        plt.clf()
        plt.close("all")

    def plot_compare_Transfer(self):
        """
        Plots and saves the comparison of all the transfer functions in a unique image
        """

        # protection from direct calls if transfer functions are not included
        if not self.transfer: return

        data1 = self.transfer_func_1
        data2 = self.transfer_func_2

        # number of transfer functions:
        num1     = data1.shape[1]-1
        num2     = data2.shape[1]-1
        # protection against different runs
        if num1!=num2:
            print( 'wrong number of transfer functions' )
            return

        # get the x values:
        xvalues_1 = data1[:,0]
        xvalues_2 = data2[:,0]
        # get the k grid:
        k_min = np.amax( [np.amin(xvalues_1), np.amin(xvalues_2)] )
        k_max = np.amin( [np.amax(xvalues_1), np.amax(xvalues_2)] )
        k_values = np.logspace(np.log10(k_min), np.log10(k_max), 1000)

        # set up the plots:
        plots_1       = CMB_plots()
        plots_2       = CMB_plots()
        plots_compa   = CMB_plots()

        plots_1.color          = self.color1
        plots_2.color          = self.color2
        plots_compa.color      = self.color_compa
        plots_compa.comparison = True
        plots_compa.axes_label_position = 'right'
        fig = plt.gcf()

        if self.transfer_labels is None:
            labels = [ 'CDM', 'baryons', 'photons', 'massless neutrinos', 'massive neutrinos',
                       'CDM+baryons+massive neutrinos', 'CDM+baryons', 'CDM+baryons+massive neutrinos+ de',
                       'Weyl potential', 'vel_Newt_cdm', 'vel_Newt_b', 'relative baryon-CDM velocity',
                      ]
        else:
            labels = self.transfer_labels

        if not len(labels) == num1:
            print( 'Not enough transfer functions for labels' )
            print( 'Labels are:' )
            print( labels )
            return

        # do the plots:
        for ind in range(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            # get the y values:
            yvalues_1    = data1[:,ind]
            yvalues_2    = data2[:,ind]
            # interpolate the two spectra:
            spline_1 = UnivariateSpline( xvalues_1, yvalues_1, s=0)
            spline_2 = UnivariateSpline( xvalues_2, yvalues_2, s=0)
            # evaluate on the l grid:
            yvalues_1 = spline_1(k_values)
            yvalues_2 = spline_2(k_values)
            # protect against zeroes:
            yvalues_1_temp = np.abs( yvalues_1 )
            yvalues_2_temp = np.abs( yvalues_2 )
            try:
                min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
            except:
                min2val_1      = cutoff
                min2val_2      = cutoff
            np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
            np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
            # computation of the percentual comparison:
            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
            # protection against values too small:
            np.place(yvalues_comp, abs(yvalues_comp)<cutoff, [cutoff])

            if not ( np.all( np.abs(yvalues_1) <= cutoff) ):
                plots_1.Transfer_plot(temp, k_values, yvalues_1)
            if not ( np.all( np.abs(yvalues_2) <= cutoff) ):
                plots_2.Transfer_plot(temp, k_values, yvalues_2)
            if not ( np.all( np.abs(yvalues_comp) <= cutoff) ):
                plots_compa.Transfer_plot(temp_comp, k_values, yvalues_comp)

            temp.set_title(labels[ind-1])

        # set the size of the image
        fig.set_size_inches( self.x_size_inc, 1.61803398875*self.x_size_inc/6.*num1 )
        # set a tight layout
        fig.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)
        # set the global title
        plt.suptitle(self.name_h1+' VS '+self.name_h2+' comparison of transfer functions', fontsize=16)
        # set the global legend
        fig.legend( handles = [plots_1.Transfer_p, plots_2.Transfer_p],
                    labels  = [self.name_h1, self.name_h2],
                    loc='lower center', ncol=3 ,fancybox=True, edgecolor='k')
        # adjust the size of the plot
        fig.subplots_adjust(top=0.95, bottom=0.05)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_transfer_comp.pdf')
        plt.clf()
        plt.close("all")
