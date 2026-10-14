import numpy             as np
import matplotlib.pyplot as plt
import matplotlib.tri    as tri       

class stability_1D:
    """ 
    class that contains the necessary tools to import and plot
    the stability space of 1 parameter EFT model
    """
    
    color = 'DarkTurquoise'
    lims  = [-0.05,1.05]
    
    def __init__(self, filename, param_name, model_name):
        """
        Class constructor:
            filename   = name of the file from which to read the data. The format should be: param  stability
            param_name = name of the parameter that is being varied
            model_name = name of the model that is being considered
        """
        # load the data
        self.data = np.loadtxt(filename)
        self.data = self.data[self.data[:,0].argsort()]
        # get the number of points
        self.numpoints = len(self.data)
        # arrange the points so that they can be digested by plot
        self.x = self.data[0:self.numpoints,0]
        self.y = self.data[0:self.numpoints,1]
        # set the param name and the model name
        self.param_name = param_name
        self.model_name = model_name

    def stab_plot(self):
        plt.plot(self.x, self.y,self.color,label='stable region')
        plt.xlabel(self.param_name, fontsize=18)
        plt.title('Stability of '+self.model_name, fontsize=14)
        plt.ylim(self.lims)
        plt.yticks([0,1])
        plt.fill_between(self.x, self.y, 0, color=self.color, alpha=0.3)
        lgd = plt.legend(loc='upper left', fancybox=True, ncol=1)
    
    
class stability_2D:
    """
    class that contains the necessary tools to import and plot
    the stability space of 2 parameters EFT model
    """
    
    color = 'DarkTurquoise'
    
    def __init__(self, filename, param_names, model_name):
        """
        Class constructor:
            filename   = name of the file from which to read the data. The format should be: param  stability
            param_name = name of the parameter that is being varied
            model_name = name of the model that is being considered
        """
        # load the data
        self.data = np.loadtxt(filename)
        # get the number of points
        self.numpoints = len(self.data)
        # arrange the points so that they can be digested by plot
        self.x = self.data[0:self.numpoints,0]
        self.y = self.data[0:self.numpoints,1]
        self.z = self.data[0:self.numpoints,2]
        self.param_names = param_names
        self.model_name  = model_name
    
    def stab_plot(self):
        x = self.x
        y = self.y
        z = self.z
        # initialize the delauney triangulation:
        triang = tri.Triangulation(x, y)
        # do the plots:
        plt.tricontourf(triang, z, colors= self.color, alpha=0.3, levels=[1., 2.])
        plt.tricontour(triang, z, colors=self.color, levels=[1., 2.])
        # title of the plot
        plt.title('Stability of '+self.model_name, fontsize=16)
        # labels
        plt.xlabel(self.param_names[0], fontsize=18)
        plt.ylabel(self.param_names[1], fontsize=18)
    
    
    
    
    
    
    
    
    
    
    
    