import numpy             as np
import matplotlib.pyplot as plt
import math

class CMB_plots:
    """
    Class that contains the methods to optimally plot the CMB power spectra
    """
    
    # general plot settings:
    color               = 'red'        # default color of the plot
    axes_label_position = 'left'       # position of the y axes and y axes label
    negative_line_style = '--'
    # comparison plot settings:
    comparison     = False             # if the plot is a comparison of spectra or just the plot
    comparison_min = 10.0**(-3)        # minimum y value of comparison plots
    comparison_max = 1.1*10.0**(+3)    # maximum y value of comparison plots
    Fsky           = 0.85              # fsky for cosmic variance

    
    def CosmicVariance(self,l):
        """ function that computes cosmic variance at a given l""" 
        return math.sqrt(2.0/((2.0*l + 1.0)*self.Fsky))
    
    def TT_plot(self, stream, xval, yval):
        """ CMB temperature power spectrum plot """   
        # do the plot:
        self.TT_p, = stream.plot( xval, yval, color = self.color )
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        # set labels
        stream.set_xlabel(r'$l$')
        stream.set_ylabel(r'$l(l+1)C_l^{TT}/ 2\pi$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            # plot cosmic variance
            ycosmicvar = np.array(map(self.CosmicVariance,xval))*100
            self.CV_p, = stream.plot(xval, ycosmicvar, color = 'k')
            self.TT_p, = stream.plot( xval, -yval, color = self.color, linestyle=self.negative_line_style )
            # set log scale
            stream.set_yscale('Log')
            # set limits and label
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta C_l^{TT}/ C_l^{TT} (\%) $')

        
    def EE_plot(self,stream, xval, yval):
        """ CMB E mode polarization power spectrum plot """
        # do the plot:
        self.EE_p, = stream.plot(xval, yval, color = self.color)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.set_yscale('Log')
        # set labels
        stream.set_xlabel(r'$l$')
        stream.set_ylabel(r'$l(l+1)C_l^{EE}/ 2\pi$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            ycosmicvar = np.array(map(self.CosmicVariance,xval))*100
            self.EE_p, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
            self.CV_p, = stream.plot(xval, ycosmicvar, color = 'k')
            stream.set_yscale('Log')
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta C_l^{EE}/ C_l^{EE} (\%) $')
                        
    def TE_plot(self,stream, xval, yval):
        """ CMB temperature E mode polarization cross correlation power spectrum plot """
        # do the plot:
        self.TE_p, = stream.plot(xval, yval, color = self.color)
        self.TE_p, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.set_yscale('Log')
        # set labels
        stream.set_xlabel(r'$l$')
        stream.set_ylabel(r'$l(l+1)C_l^{TE}/ 2\pi$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            ycosmicvar = np.array(map(self.CosmicVariance,xval))*100
            self.CV_p, = stream.plot(xval, ycosmicvar, color = 'k')
            stream.set_yscale('Log')
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta C_l^{TE}/ C_l^{TE} (\%) $')
            
    def BB_plot(self,stream, xval, yval):
        """ CMB B mode polarization power spectrum plot """
        # do the plot:
        self.BB_p, = stream.plot(xval, yval, color = self.color)
        self.BB_p, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.set_yscale('Log')
        # set labels
        stream.set_xlabel(r'$l$')
        stream.set_ylabel(r'$l(l+1)C_l^{BB}/ 2\pi$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            ycosmicvar = np.array(map(self.CosmicVariance,xval))*100
            self.CV_p, = stream.plot(xval, ycosmicvar, color = 'k')
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta C_l^{BB}/ C_l^{BB} (\%) $')    
                                
    def Phi_plot(self,stream, xval, yval):
        """ CMB lensing power spectrum plot """
        # do the plot:
        self.Phi_p, = stream.plot(xval, yval, color = self.color)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        # set labels
        stream.set_xlabel(r'$l$')
        stream.set_ylabel(r'$l^4 C_l^{\Phi\Phi}$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            ycosmicvar = np.array(map(self.CosmicVariance,xval))*100
            self.Phi_p, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
            self.CV_p, = stream.plot(xval, ycosmicvar, color = 'k')
            stream.set_yscale('Log')
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta C_l^{\Phi\Phi}/ C_l^{\Phi\Phi} (\%) $') 
                        
    def PhiT_plot(self,stream, xval, yval):
        """ CMB lensing and temperature cross correlation power spectrum plot """
        # do the plot:
        self.PhiT_p, = stream.plot(xval, yval, color = self.color)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        stream.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        # set axes scales
        stream.set_xscale('Log')
        # set labels
        stream.set_xlabel(r'$l$')
        stream.set_ylabel(r'$l^3 C_l^{\Phi T}$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            ycosmicvar = np.array(map(self.CosmicVariance,xval))*100
            self.PhiT_p, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
            self.CV_p, = stream.plot(xval, ycosmicvar, color = 'k')
            stream.set_yscale('Log')
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta C_l^{\Phi T}/ C_l^{\Phi T} (\%) $')  
    
    def Generic_Cl(self,stream, xval, yval):
        """ Generic spectrum plotter (in l-space) """
        # take the abs value of y-val
        yval = np.array(yval)
        # do the plot:
        self.Generic_Cl_plot, = stream.plot(xval, yval, color = self.color)
        self.Generic_Cl_plot, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.set_yscale('Log')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            ycosmicvar = np.array(map(self.CosmicVariance,xval))*100
            self.CV_p, = stream.plot(xval, ycosmicvar, color = 'k')
            stream.set_ylim(self.comparison_min, self.comparison_max)
  
                            
    def Matter_plot(self,stream, xval, yval):
        """ Matter power spectrum plot """
        # do the plot:
        self.Matter_p, = stream.plot(xval, yval, color = self.color)
        self.Matter_p, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.set_yscale('Log')
        # set labels
        stream.set_xlabel(r'$k$')
        stream.set_ylabel(r'$P(k)$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:    
        if self.comparison:
            stream.set_yscale('Log')
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta P(k)/ P(k) (\%) $')

    def Transfer_plot(self,stream, xval, yval):
        """ Transfer functions plot """
        # do the plot:
        self.Transfer_p, = stream.plot(xval, yval, color = self.color)
        self.Transfer_p, = stream.plot(xval, -yval, color = self.color, linestyle=self.negative_line_style)
        # set x axes boundary:
        stream.set_xlim(np.amin(xval),np.amax(xval))
        # set axes scales
        stream.set_xscale('Log')
        stream.set_yscale('Log')
        # set labels
        stream.set_ylabel(r'$T(k)$')
        # set the position of axes and label:
        stream.yaxis.set_label_position(self.axes_label_position)
        if self.axes_label_position == 'right':
            stream.yaxis.tick_right()
        # setup if comparison:
        if self.comparison:
            stream.set_yscale('Log')
            stream.set_ylim(self.comparison_min, self.comparison_max)
            stream.set_ylabel(r'$\Delta T(k)/ T(k) (\%) $')







 