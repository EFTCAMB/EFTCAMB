import numpy             as np
import matplotlib.pyplot as plt
import math

from CMB_plots import CMB_plots


class CAMB_results_plot:
    """ 
    class that contains the necessary tools to plot the results of a CAMB run
    """
    
    # general plot settings:
    color      = 'red'                 # default color of the plot
    x_size_inc = 8.30                  # x size of the final plots, in inches
    y_size_inc = 8.30                  # y size of the final plots, in inches   
    
    def __init__(self, root, outpath, tensor=False, lensing=False, transfer=False, name=''):
        """
        Class constructor:
            root     = name of the CAMB run
            outpath  = path to the output folder
            tensor   = specifies wether the results contains the tensor Cls
            lensing  = specifies wether the results contains lensing
            transfer = specifies wether the results contains transfer functions
        """
        
        # store the constructor options:
        self.root     = root
        self.outpath  = outpath
        self.tensor   = tensor
        self.lensing  = lensing
        self.transfer = transfer
        
        # extract the name of the model from the root:
        self.name   = root.split("/")[len(root.split("/"))-1]
        
        # set the human readable name
        if name=='':
            self.name_h = ''.join(i for i in self.name.replace('_',' ') if not i.isdigit())
        else:
            self.name_h = name
        
        # load the data:
        self.scalCls           = np.loadtxt(root+'_scalCls.dat')
        self.scalCovCls        = np.loadtxt(root+'_scalCovCls.dat')
        
        if self.lensing:
            self.lensedCls         = np.loadtxt(root+'_lensedCls.dat')
            self.lensedtotCls      = np.loadtxt(root+'_lensedtotCls.dat')
            self.lenspotentialCls  = np.loadtxt(root+'_lenspotentialCls.dat')
        
        if self.transfer:
            self.matterpower       = np.loadtxt(root+'_matterpower.dat')
            self.transfer_func     = np.loadtxt(root+'_transfer_out.dat')
        
        if self.tensor:
            self.tensCls           = np.loadtxt(root+'_tensCls.dat')
            self.totCls            = np.loadtxt(root+'_totCls.dat')

        
    def plot_scalCls(self):
        """
        Plots and saves all the scalar Cls in a unique image
        """
        
        # number of Cls:
        num     = self.scalCls.shape[1]-1
        if self.transfer:
            num += 1
        
        # l values:
        xvalues = self.scalCls[:,0]
        # set up the plots:
        plots       = CMB_plots()
        plots.color = self.color
        
        index = 0
        fig = plt.gcf()
        
        for ind in xrange(1,num+1):
               
            index2 = int((ind-1)/2)
            
            temp = plt.subplot2grid((int(math.ceil(num/2.0)),2),(index2, index))
            
            if ind==1:
                yvalues = self.scalCls[:,ind]
                plots.TT_plot(temp, xvalues, yvalues)
                temp.set_title('TT power spectrum')
            elif ind == 2:
                yvalues = self.scalCls[:,ind]
                plots.EE_plot(temp, xvalues, yvalues)
                temp.set_title('EE power spectrum')
            elif ind == 3:
                yvalues = self.scalCls[:,ind]
                plots.TE_plot(temp, xvalues, yvalues)
                temp.set_title('TE power spectrum')
            elif ind == 4 and self.lensing:
                yvalues = self.scalCls[:,ind]
                plots.Phi_plot(temp, xvalues, yvalues)
                temp.set_title('$\phi$ power spectrum')
            elif ind == 5 and self.lensing:
                yvalues = self.scalCls[:,ind]
                plots.PhiT_plot(temp, xvalues, yvalues)
                temp.set_title('$\phi$T power spectrum')
            elif ind == num and self.transfer:
                xvalues = self.matterpower[:,0]
                yvalues = self.matterpower[:,1]
                plots.Matter_plot(temp, xvalues, yvalues)
                temp.set_title('Matter power spectrum')
            else:
                yvalues = self.scalCls[:,ind]
                plots.Generic_Cl(temp, xvalues, yvalues)
            
            if index==0: index = 1
            elif index==1: index = 0

        # global title of the plot
        plt.suptitle(self.name_h+' scalar Cls', fontsize=16)
        # define the image size and layout
        fig.set_size_inches(self.x_size_inc, self.y_size_inc)
        fig.tight_layout()
        fig.subplots_adjust(top=0.90) # adjust for the fact that tight_layout does not consider suptitle 
        # save the figure and close
        plt.savefig(self.outpath+self.name+'_scalCls.png')
        plt.clf()
    
    
    def plot_scalCovCls(self):
        """
        Plots and saves all the scalar covariance Cls in a unique image
        """
        # number of Cls:
        num     = int(math.sqrt(self.scalCovCls.shape[1]-1))
        
        # l values:
        xvalues = self.scalCovCls[:,0]
        # set up the plots:
        plots       = CMB_plots()
        plots.color = self.color
        
        fig = plt.gcf()
        
        # setup a dictionary with the names of the Cls
        dict = { 1: 'T', 2: 'E', 3: '$\phi$'}
        for i in xrange(4, num+1):
            dict[i] = 'W'+str(i)
            
        for ind in xrange(1, num+1):
            for ind2 in xrange(1, ind+1):
                
                temp = plt.subplot2grid((num, num),(ind-1, ind2-1))
                col = ind + num*(ind2-1)
                yvalues = self.scalCovCls[:,col]
                
                temp.text( 0.15, 0.15,  dict[ind]+dict[ind2] , ha='center', va='center', transform=temp.transAxes)
                
                plots.Generic_Cl(temp, xvalues, yvalues)
                
                if ind == num:
                    temp.set_xlabel(r'$l$')
                
        # global title of the plot
        plt.suptitle(self.name_h+' scalar Cov Cls', fontsize=16)
        # define the image size and layout. We need a bigger image for the covariance.
        fig.set_size_inches(1.61803398875*self.x_size_inc, self.x_size_inc)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0)
        fig.subplots_adjust(top=0.90) # adjust for the fact that tight_layout does not consider suptitle 
        # save the figure and close
        plt.savefig(self.outpath+self.name+'_scalCovCls.png')
        plt.clf()
    
    def plot_lensedCls(self):
        """
        Plots and saves all the lensed Cls in a unique image
        """
        
        # protection from direct calls if lensing is not included
        if not self.lensing: return
        
        # number of Cls:
        num     = self.lensedCls.shape[1]-1
        
        # l values:
        xvalues = self.lensedCls[:,0]
        # set up the plots:
        plots       = CMB_plots()
        plots.color = self.color
        
        index = 0
        fig = plt.gcf()
        
        for ind in xrange(1,num+1):

            yvalues = self.lensedCls[:,ind]   
            index2 = int((ind-1)/2)
            
            temp = plt.subplot2grid((int(math.ceil(num/2.0)),2),(index2, index))
            
            if ind==1:
                plots.TT_plot(temp, xvalues, yvalues)
                temp.set_title('TT power spectrum')
            elif ind == 2:
                plots.EE_plot(temp, xvalues, yvalues)
                temp.set_title('EE power spectrum')
            elif ind == 3:
                plots.BB_plot(temp, xvalues, yvalues)
                temp.set_title('BB power spectrum')
            elif ind == 4:
                plots.TE_plot(temp, xvalues, yvalues)
                temp.set_title('TE power spectrum')
            else:
                temp.plot(xvalues,yvalues, color = self.color)
            
            if index==0: index = 1
            elif index==1: index = 0
        
        # global title of the plot
        plt.suptitle(self.name_h+' lensed Cls', fontsize=16)
        # define the image size and layout
        fig.set_size_inches(self.x_size_inc, self.y_size_inc)
        fig.tight_layout()
        fig.subplots_adjust(top=0.90) # adjust for the fact that tight_layout does not consider suptitle
        # save the figure and close
        plt.savefig(self.outpath+self.name+'_lensCls.png')  
        plt.clf()
 
 
    def plot_tensCls(self):
        """
        Plots and saves all the tensor Cls in a unique image
        """
        
        # protection from direct calls if tensors is not included
        if not self.tensor: return
        
        # number of Cls:
        num     = self.tensCls.shape[1]-1
        
        # l values:
        xvalues = self.tensCls[:,0]
        # set up the plots:
        plots       = CMB_plots()
        plots.color = self.color
        
        index = 0
        fig = plt.gcf()
        
        for ind in xrange(1,num+1):
            
            yvalues = self.tensCls[:,ind]
            index2 = int((ind-1)/2)
            
            temp = plt.subplot2grid((int(math.ceil(num/2.0)),2),(index2, index))
            
            if ind==1:
                plots.TT_plot(temp, xvalues, yvalues)
                temp.set_yscale('Log')
                temp.set_title('TT power spectrum')
            elif ind == 2:
                plots.EE_plot(temp, xvalues, yvalues)
                temp.set_title('EE power spectrum')
            elif ind == 3:
                plots.BB_plot(temp, xvalues, yvalues)
                temp.set_title('BB power spectrum')
            elif ind == 4:
                plots.TE_plot(temp, xvalues, yvalues)
                temp.set_title('TE power spectrum')
            else:
                temp.plot(xvalues,yvalues, color = self.color)
            
            if index==0: index = 1
            elif index==1: index = 0

        # global title of the plot
        plt.suptitle(self.name_h+' tensor Cls', fontsize=16)
        # define the image size and layout
        fig.set_size_inches(self.x_size_inc, self.y_size_inc)
        fig.tight_layout()
        fig.subplots_adjust(top=0.90)
        # save the figure and close
        plt.savefig(self.outpath+self.name+'_tensCls.png')
        plt.clf()


    def plot_totalCls(self):
        """
        Plots and saves all the total (scalar + tensor) Cls in a unique image
        If lensing is included the results are lensed Cls.
        """
        
        # protection from direct calls if tensors is not included
        if not self.tensor: return
        
        # number of Cls:
        num     = self.totCls.shape[1]-1
        if self.lensing:
            num = self.lensedtotCls.shape[1]-1
        
        # l values:
        xvalues = self.totCls[:,0]
        if self.lensing:
            xvalues = self.lensedtotCls[:,0]
        
        # set up the plots:
        plots       = CMB_plots()
        plots.color = self.color
        
        index = 0
        fig = plt.gcf()
        
        for ind in xrange(1,num+1):
            
            yvalues = self.totCls[:,ind]
            if self.lensing:
                yvalues = self.lensedtotCls[:,ind]
            
            index2 = int((ind-1)/2)
            
            temp = plt.subplot2grid((int(math.ceil(num/2.0)),2),(index2, index))
            
            if ind==1:
                plots.TT_plot(temp, xvalues, yvalues)
                temp.set_title('TT power spectrum')
            elif ind == 2:
                plots.EE_plot(temp, xvalues, yvalues)
                temp.set_title('EE power spectrum')
            elif ind == 3:
                plots.BB_plot(temp, xvalues, yvalues)
                temp.set_title('BB power spectrum')
            elif ind == 4:
                plots.TE_plot(temp, xvalues, yvalues)
                temp.set_title('TE power spectrum')
            else:
                temp.plot(xvalues,yvalues, color = self.color)
            
            if index==0: index = 1
            elif index==1: index = 0

        # global title of the plot
        if self.lensing:
            plt.suptitle(self.name_h+' total lensed Cls', fontsize=16)
        else:
            plt.suptitle(self.name_h+' total Cls', fontsize=16)
        # define the image size and layout
        fig.set_size_inches(self.x_size_inc, self.y_size_inc)
        fig.tight_layout()
        fig.subplots_adjust(top=0.90) # adjust for the fact that tight_layout does not consider suptitle
        # save the figure and close
        plt.savefig(self.outpath+self.name+'_totCls.png')
        plt.clf()


    def plot_Transfer(self):
        """
        Plots and saves all the transfer functions in a unique image
        """
        
        # protection from direct calls if transfer functions are not included
        if not self.transfer: return
    
        # number of transfer functions:
        num     = self.transfer_func.shape[1]-1
    
        # k values:
        xvalues = self.transfer_func[:,0]

        # set up the plots:
        plots       = CMB_plots()
        plots.color = self.color

        index = 0
        fig = plt.gcf()
        
        labels = [ 'CDM', 'baryons', 'photons', 'massless neutrinos', 'massive neutrinos',
                   'CDM+baryons+massive neutrinos', 'CDM+baryons', 'CDM+baryons+massive neutrinos+ de',
                   'The Weyl potential', 'vel_Newt_cdm', 'vel_Newt_b', 'relative baryon-CDM velocity'
                  ]
        
        for ind in xrange(1,num+1):

            yvalues = self.transfer_func[:,ind]
            index2 = int((ind-1)/2)
            
            temp = plt.subplot2grid((int(math.ceil(num/2.0)),2),(index2, index))
            
            plots.Transfer_plot(temp, xvalues, yvalues)
            temp.set_title(labels[ind-1])
            
            if index==0: index = 1
            elif index==1: index = 0
            
            if ind == num or ind == num-1:
                temp.set_xlabel(r'$k$')
                
        plt.suptitle(self.name_h+' transfer functions', fontsize=16)
        # define the image size and layout
        fig.set_size_inches(self.x_size_inc, 1.61803398875*self.x_size_inc)
        fig.tight_layout()
        fig.subplots_adjust(top=0.94) # adjust for the fact that tight_layout does not consider suptitle
        # save the figure and close
        plt.savefig(self.outpath+self.name+'_transfer.png')
        plt.clf()  
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        