import numpy             as np
import matplotlib.pyplot as plt
import math

from CMB_plots import CMB_plots


class CAMB_results_compare_plot:
    """
    class that contains the necessary tools to plot the comparison of two CAMB run
    """

    # general plot settings:
    color1      = 'red'                # default color of the line of the first model
    color2      = 'blue'               # default color of the line of the second model
    color_compa = 'green'              # default color of the line of the difference
    x_size_inc  = 8.30                 # x size of the final plots, in inches
    y_size_inc  = 11.7                 # y size of the final plots, in inches

    def __init__(self, root1, root2, outpath,
                 tensor=False, lensing=False, transfer=False,
                 name1='', name2=''):
        """
        Class constructor:
            root1     = name of the first CAMB run
            root2     = name of the second CAMB run
            outpath   = path to the output folder
            tensor    = (optional) specifies wether the results contains the tensor Cls
            lensing   = (optional) specifies wether the results contains lensing
            transfer  = (optional) specifies wether the results contains transfer functions
            name1     = (optional) specifies the name of the first model. Used for the legend.
            name2     = (optional) specifies the name of the second model. Used for the legend.
        """

        # store the constructor options:
        self.root1    = root1
        self.root2    = root2
        self.outpath  = outpath
        self.tensor   = tensor
        self.lensing  = lensing
        self.transfer = transfer

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
        self.scalCovCls_1        = np.loadtxt(root1+'_scalCovCls.dat')

        self.scalCls_2           = np.loadtxt(root2+'_scalCls.dat')
        self.scalCovCls_2        = np.loadtxt(root2+'_scalCovCls.dat')

        if self.lensing:
            self.lensedCls_1         = np.loadtxt(root1+'_lensedCls.dat')
            self.lenspotentialCls_1  = np.loadtxt(root1+'_lenspotentialCls.dat')

            self.lensedCls_2         = np.loadtxt(root2+'_lensedCls.dat')
            self.lenspotentialCls_2  = np.loadtxt(root2+'_lenspotentialCls.dat')

        if self.transfer:
            self.matterpower_1       = np.loadtxt(root1+'_matterpower.dat')
            self.transfer_func_1     = np.loadtxt(root1+'_transfer_out.dat')

            self.matterpower_2       = np.loadtxt(root2+'_matterpower.dat')
            self.transfer_func_2     = np.loadtxt(root2+'_transfer_out.dat')

        if self.tensor:
            self.tensCls_1           = np.loadtxt(root1+'_tensCls.dat')
            self.totCls_1            = np.loadtxt(root1+'_totCls.dat')

            self.tensCls_2           = np.loadtxt(root2+'_tensCls.dat')
            self.totCls_2            = np.loadtxt(root2+'_totCls.dat')

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
            print 'wrong number of Cls'
            return

        if len(self.scalCls_1[:,0])!=len(self.scalCls_2[:,0]):
            print 'different lmax'
            return

        if len(self.matterpower_1[:,0])!=len(self.matterpower_2[:,0]):
            self.transfer = False

        # add the matter power spectrum if required:
        if self.transfer: num1 += 1
        # values of l
        xvalues = self.scalCls_1[:,0]

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
        for ind in xrange(1,num1+1):
            # distribute the plots in the figure:
            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            if ind == 1: # TT power spectrum:
                yvalues_1    = self.scalCls_1[:,ind]
                yvalues_2    = self.scalCls_2[:,ind]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.TT_plot(temp, xvalues, yvalues_1)
                plots_2.TT_plot(temp, xvalues, yvalues_2)
                plots_compa.TT_plot(temp_comp, xvalues, yvalues_comp)
                temp_comp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2: # EE power spectrum:
                yvalues_1    = self.scalCls_1[:,ind]
                yvalues_2    = self.scalCls_2[:,ind]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.EE_plot(temp, xvalues, yvalues_1)
                plots_2.EE_plot(temp, xvalues, yvalues_2)
                plots_compa.EE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3: # TE power spectrum:
                yvalues_1    = self.scalCls_1[:,ind]
                yvalues_2    = self.scalCls_2[:,ind]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.TE_plot(temp, xvalues, yvalues_1)
                plots_2.TE_plot(temp, xvalues, yvalues_2)
                plots_compa.TE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('TE power spectrum')

            elif ind == 4 and self.lensing: # CMB lensing power spectrum:
                yvalues_1    = self.scalCls_1[:,ind]
                yvalues_2    = self.scalCls_2[:,ind]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.Phi_plot(temp, xvalues, yvalues_1)
                plots_2.Phi_plot(temp, xvalues, yvalues_2)
                plots_compa.Phi_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('$\phi$ power spectrum')

            elif ind == 5 and self.lensing: # CMB lensing - Temperature power spectrum:
                yvalues_1    = self.scalCls_1[:,ind]
                yvalues_2    = self.scalCls_2[:,ind]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.PhiT_plot(temp, xvalues, yvalues_1)
                plots_2.PhiT_plot(temp, xvalues, yvalues_2)
                plots_compa.PhiT_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('$\phi$T power spectrum')

            elif ind == num1 and self.transfer: # matter power spectrum:
                xvalues      = self.matterpower_2[:,0]
                yvalues_1    = self.matterpower_1[:,1]
                yvalues_2    = self.matterpower_2[:,1]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.Matter_plot(temp, xvalues, yvalues_1)
                plots_2.Matter_plot(temp, xvalues, yvalues_2)
                plots_compa.Matter_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('Matter power spectrum')

            else: # generic Cl comparison:
                yvalues_1    = self.scalCls_1[:,ind]
                yvalues_2    = self.scalCls_2[:,ind]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.Generic_Cl(temp, xvalues, yvalues_1)
                plots_2.Generic_Cl(temp, xvalues, yvalues_2)
                plots_compa.Generic_Cl(temp_comp, xvalues, yvalues_comp)

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
                    loc='lower center', ncol=3 ,fancybox=True)
        # adjust the size of the plot
        fig.subplots_adjust(top=0.92, bottom=0.08)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_scalCls_comp.pdf')
        plt.clf()
        plt.close("all")

    def plot_compare_scalCovCls(self):
        """
        Plots and saves the comparison of all the scalar Cov Cls in a unique image
        """

        # number of Cls:
        num1     = self.scalCovCls_1.shape[1]-1
        num2     = self.scalCovCls_2.shape[1]-1
        # protection against different runs
        if num1!=num2:
            print 'wrong number of Cls'
            return

        if len(self.scalCovCls_1[:,0])!=len(self.scalCovCls_2[:,0]):
            print 'different lmax'
            return

        # size of the Cl Cov matrix:
        num1     = int(math.sqrt(num1))
        num_tot  = sum(xrange(1,num1+1))

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
        for i in xrange(4, num1+1):
            dict[i] = 'W'+str(i)

        # values of the multipoles:
        xvalues = self.scalCovCls_1[:,0]
        # other stuff:
        ind_tot = 0
        # do the plots:
        for ind in xrange(1,num1+1):
            for ind2 in xrange(1, ind+1):

                ind_tot += 1
                # place the plots:
                temp      = plt.subplot2grid((num_tot,2), (ind_tot-1, 0))
                temp_comp = plt.subplot2grid((num_tot,2), (ind_tot-1, 1))
                # values of the Cls:
                col          = ind + num1*(ind2-1)
                yvalues_1    = self.scalCovCls_1[:,col]
                yvalues_2    = self.scalCovCls_2[:,col]
                # protection against values equal to zero:
                yvalues_temp = np.array(map(abs, yvalues_1))
                min2val      = np.min(yvalues_temp[np.nonzero(yvalues_temp)])
                np.place(yvalues_1, yvalues_1 == 0.0, min2val)
                # computation of the percentual comparison:
                yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
                # protection against values too small:
                np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])
                # make the plots:
                plots_1.Generic_Cl(temp, xvalues, yvalues_1)
                plots_2.Generic_Cl(temp, xvalues, yvalues_2)
                plots_compa.Generic_Cl(temp_comp, xvalues, yvalues_comp)

                temp.set_title(dict[ind2]+dict[ind]+' power spectrum')

        # set the size of the image
        fig.set_size_inches( self.x_size_inc, self.y_size_inc/5.0*num_tot)
        # set a tight layout
        fig.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)

        # set the global legend
        fig.legend( handles = [plots_1.Generic_Cl_plot, plots_2.Generic_Cl_plot, plots_compa.CV_p],
                    labels  = [self.name_h1, self.name_h2, 'Cosmic variance'],
                    loc='lower center', ncol=3 ,fancybox=True)

        # set the global title
        plt.suptitle(self.name_h1+' VS '+self.name_h2+
                     ' comparison of scalar Cov Cls', fontsize=16)

        # adjust the size of the plot
        fig.subplots_adjust(top=0.92, bottom=0.08)
        #        fig.subplots_adjust(top=0.96, bottom=0.01)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_scalCovCls_comp.pdf')
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
            print 'wrong number of Cls'
            return

        if len(self.lensedCls_1[:,0])!=len(self.lensedCls_2[:,0]):
            print 'different lmax'
            return

        xvalues = self.lensedCls_1[:,0]

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
        for ind in xrange(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            yvalues_1    = self.lensedCls_1[:,ind]
            yvalues_2    = self.lensedCls_2[:,ind]

            min2val      = np.min(yvalues_1[np.nonzero(yvalues_1)])
            np.place(yvalues_1, yvalues_1 == 0.0, min2val)

            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
            # Protection against all zero: put to machine precision the values that are zero
            np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])

            if ind == 1:
                plots_1.TT_plot(temp, xvalues, yvalues_1)
                plots_2.TT_plot(temp, xvalues, yvalues_2)
                plots_compa.TT_plot(temp_comp, xvalues, yvalues_comp)
                temp_comp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2:
                plots_1.EE_plot(temp, xvalues, yvalues_1)
                plots_2.EE_plot(temp, xvalues, yvalues_2)
                plots_compa.EE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3:
                plots_1.BB_plot(temp, xvalues, yvalues_1)
                plots_2.BB_plot(temp, xvalues, yvalues_2)
                plots_compa.BB_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('BB power spectrum')

            elif ind == 4:
                plots_1.TE_plot(temp, xvalues, yvalues_1)
                plots_2.TE_plot(temp, xvalues, yvalues_2)
                plots_compa.TE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('TE power spectrum')

            else:
                plots_1.TT_plot(temp, xvalues, yvalues_1)
                plots_2.TT_plot(temp, xvalues, yvalues_2)
                plots_compa.TT_plot(temp_comp, xvalues, yvalues_comp)

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
                    loc='lower center', ncol=3 ,fancybox=True)
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
            print 'wrong number of Cls'
            return

        if len(self.tensCls_1[:,0])!=len(self.tensCls_2[:,0]):
            print 'different lmax'
            return

        xvalues = self.tensCls_1[:,0]

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
        for ind in xrange(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            yvalues_1    = self.tensCls_1[:,ind]
            yvalues_2    = self.tensCls_2[:,ind]

            min2val      = np.min(yvalues_1[np.nonzero(yvalues_1)])
            np.place(yvalues_1, yvalues_1 == 0.0, min2val)

            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
            # Protection against all zero: put to machine precision the values that are zero
            np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])

            if ind == 1:
                plots_1.TT_plot(temp, xvalues, yvalues_1)
                plots_2.TT_plot(temp, xvalues, yvalues_2)
                plots_compa.TT_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2:
                plots_1.EE_plot(temp, xvalues, yvalues_1)
                plots_2.EE_plot(temp, xvalues, yvalues_2)
                plots_compa.EE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3:
                plots_1.BB_plot(temp, xvalues, yvalues_1)
                plots_2.BB_plot(temp, xvalues, yvalues_2)
                plots_compa.BB_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('BB power spectrum')

            elif ind == 4:
                plots_1.TE_plot(temp, xvalues, yvalues_1)
                plots_2.TE_plot(temp, xvalues, yvalues_2)
                plots_compa.TE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('TE power spectrum')

            else:
                plots_1.TT_plot(temp, xvalues, yvalues_1)
                plots_2.TT_plot(temp, xvalues, yvalues_2)
                plots_compa.TT_plot(temp_comp, xvalues, yvalues_comp)


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
                    loc='lower center', ncol=3 ,fancybox=True)
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
            print 'wrong number of Cls'
            return

        if len(data1[:,0])!=len(data2[:,0]):
            print 'different lmax'
            return

        xvalues = data1[:,0]

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
        for ind in xrange(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            yvalues_1    = data1[:,ind]
            yvalues_2    = data2[:,ind]

            min2val      = np.min(yvalues_1[np.nonzero(yvalues_1)])
            np.place(yvalues_1, yvalues_1 == 0.0, min2val)

            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
            # Protection against all zero: put to machine precision the values that are zero
            np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])

            if ind == 1:
                plots_1.TT_plot(temp, xvalues, yvalues_1)
                plots_2.TT_plot(temp, xvalues, yvalues_2)
                plots_compa.TT_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_yscale('Log')
                temp.set_title('TT power spectrum')

            elif ind == 2:
                plots_1.EE_plot(temp, xvalues, yvalues_1)
                plots_2.EE_plot(temp, xvalues, yvalues_2)
                plots_compa.EE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('EE power spectrum')

            elif ind == 3:
                plots_1.BB_plot(temp, xvalues, yvalues_1)
                plots_2.BB_plot(temp, xvalues, yvalues_2)
                plots_compa.BB_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('BB power spectrum')

            elif ind == 4:
                plots_1.TE_plot(temp, xvalues, yvalues_1)
                plots_2.TE_plot(temp, xvalues, yvalues_2)
                plots_compa.TE_plot(temp_comp, xvalues, yvalues_comp)
                temp.set_title('TE power spectrum')

            else:
                plots_1.TT_plot(temp, xvalues, yvalues_1)
                plots_2.TT_plot(temp, xvalues, yvalues_2)
                plots_compa.TT_plot(temp_comp, xvalues, yvalues_comp)


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
                    loc='lower center', ncol=3 ,fancybox=True)
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
            print 'wrong number of transfer functions'
            return
        if len(data1[:,0])!=len(data2[:,0]):
            print 'Different values of k'
            return

        xvalues = data1[:,0]

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

        labels = [ 'CDM', 'baryons', 'photons', 'massless neutrinos', 'massive neutrinos',
                   'CDM+baryons+massive neutrinos', 'CDM+baryons', 'CDM+baryons+massive neutrinos+ de',
                   'The Weyl potential', 'vel_Newt_cdm', 'vel_Newt_b', 'relative baryon-CDM velocity'
                  ]


        if not len(labels) == num1:
            print 'Not enough transfer functions'
            return

        # do the plots:
        for ind in xrange(1,num1+1):

            temp      = plt.subplot2grid((num1,2), (ind-1, 0))
            temp_comp = plt.subplot2grid((num1,2), (ind-1, 1))

            yvalues_1    = data1[:,ind]
            yvalues_2    = data2[:,ind]

            min2val      = np.min(yvalues_1[np.nonzero(yvalues_1)])
            np.place(yvalues_1, yvalues_1 == 0.0, min2val)

            yvalues_comp = (yvalues_1 - yvalues_2)/abs(yvalues_1)*100
            # Protection against all zero: put to machine precision the values that are zero
            np.place(yvalues_comp, abs(yvalues_comp)<10.0**(-16), [10.0**(-16)])

            plots_1.Transfer_plot(temp, xvalues, yvalues_1)
            plots_2.Transfer_plot(temp, xvalues, yvalues_2)
            plots_compa.Transfer_plot(temp_comp, xvalues, yvalues_comp)

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
                    loc='lower center', ncol=3 ,fancybox=True)
        # adjust the size of the plot
        fig.subplots_adjust(top=0.95, bottom=0.05)
        # save the result and close
        plt.savefig(self.outpath+self.name1+'_'+self.name2+'_transfer_comp.pdf')
        plt.clf()
        plt.close("all")
