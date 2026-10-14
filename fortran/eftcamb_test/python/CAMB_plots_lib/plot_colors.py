class _color_scheme_1D:
    """
    Class that contains several color schemes optimized for 1D plots.
    
    A color scheme is a python list consisting of the RGB codes of 6 colors that perform well
    when used on a plot involving lines.
    """
    
    scheme_1 = [(203.0/255.0, 15.0/255.0, 40.0/255.0),
                (255.0/255.0, 165.0/255.0, 0.0),
                (0.0, 0.0, 0.0),
                (42.0/255.0, 46.0/255.0, 139.0/255.0),
                (0.0/255.0, 153.0/255.0, 204.0/255.0),
                (0.0/255.0, 221.0/255.0, 52.0/255.0)
                ]    
    
    # default color scheme
    default = scheme_1
    
    # call method that returns the default color scheme
    def __call__(self):
        return self.default
    
""" ************************************************************************************************ """

class _color_scheme_2D:
    """
    Class that contains several color schemes optimized for 2D plots.
    
    A color scheme is a python list consisting of the RGB codes of 6 colors that perform well
    when used on a plot involving 2D sets.
    """
    
    scheme_1 = [(203.0/255.0, 15.0/255.0, 40.0/255.0),
                (255.0/255.0, 165.0/255.0, 0.0),
                (0.0, 0.0, 0.0),
                (42.0/255.0, 46.0/255.0, 139.0/255.0),
                (0.0/255.0, 153.0/255.0, 204.0/255.0),
                (0.0/255.0, 221.0/255.0, 52.0/255.0)
                ]    
    
    # default color scheme
    default = scheme_1
    
    # call method that returns the default color scheme
    def __call__(self):
        return self.default    
    
""" ************************************************************************************************ """
    
class plot_colors:
    """
    Class that contains the interface between the two color classes and the rest
    """
    
    color_scheme_1D = _color_scheme_1D()
    color_scheme_2D = _color_scheme_2D()
        
    
    
    
    