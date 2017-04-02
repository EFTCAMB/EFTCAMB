#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2017 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

"""

Simple Python code to plot debug quantities.

The ouput will be a plot with the specified file.

Invoking the help option ``debug_plotter.py -h`` will result in::

    usage: debug_plotter.py [-h] [-o OUTROOT] [-x X_COLUMN]
                            [-y Y_COLUMNS [Y_COLUMNS ...]] [-xl] [-yl] [-l] [-v]
                            [-q]
                            files [files ...]
    
    debug plotter
    
    positional arguments:
      files                 a list of files with Fisher matrices
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTROOT, --outroot OUTROOT
                            path, name and format of the output file. If this
                            option is not present the plotter will try to display
                            to screen.
      -x X_COLUMN, --x X_COLUMN
                            decides which column of the file to use as the x
                            coordinate for the plot.
      -y Y_COLUMNS [Y_COLUMNS ...], --y Y_COLUMNS [Y_COLUMNS ...]
                            decides which columns of the file to use as the y
                            coordinate for the plot.
      -xl, --xlog           decides wether the x scale is a log scale. If log plotting will plot the absolute value.
      -yl, --ylog           decides wether the y scale is a log scale. If log plotting will plot the absolute value.
      -l, --latex           decides wether text is latex rendered or not
      -v, --version         show program's version number and exit
      -q, --quiet           decides wether something gets printed to screen or not

Developed by Marco Raveri (mraveri@uchicago.edu) for the EFTCAMB code.

"""

# ***************************************************************************************

__version__ = '1.0' #: version of the application

# ***************************************************************************************

""" Hard coded options """

# ***************************************************************************************

# import first dependencies:
import matplotlib
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines    as mlines
import numpy as np
import argparse
import math
import sys
import os

# get the path of the application and the CosmicFish library:
here = os.path.dirname(os.path.abspath(__file__))

# ***************************************************************************************

# protection against importing:
if __name__ == "__main__":
    
    # parse command line arguments:
    parser = argparse.ArgumentParser(description='debug plotter')
    # parse file names:
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                         help='a list of files with Fisher matrices')
    # parse the output root:
    parser.add_argument('-o','--outroot', dest='outroot', type=str,
                         help='path, name and format of the output file. If this option is not present the plotter will try to display to screen.')
    # reference column:
    parser.add_argument('-x','--x', dest='x_column', type=int,
                        help='decides which column of the file to use as the x coordinate for the plot.')
    # select y columns:
    parser.add_argument('-y','--y', dest='y_columns', type=int, nargs='+',
                        help='decides which columns of the file to use as the y coordinate for the plot.')
    # x scale log option:
    parser.add_argument('-xl','--xlog', dest='xlog', action='store_true', 
                        help='decides wether the x scale is a log scale. If log plotting will plot the absolute value.')
    # y scale log option:
    parser.add_argument('-yl','--ylog', dest='ylog', action='store_true', 
                        help='decides wether the y scale is a log scale. If log plotting will plot the absolute value.')
    # latex text option:
    parser.add_argument('-l','--latex', dest='latex', action='store_true', 
                        help='decides wether text is latex rendered or not')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
            
    # process input arguments:
    files          = args.files
    number_fish    = len(files)
    if args.outroot is not None:
        outroot    = args.outroot
    else:
        outroot    = None
    
    # latex rendering of text if required:   
    if args.latex is not None:
        if args.latex:
            plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
            plt.rc('text', usetex=True)
    
    # start the plot:
    fig = plt.figure()
    
    # cycle over files:
    for file in files:
        # import the data:
        data  = np.loadtxt( file )
        # get number of columns:
        n_col = data.shape[1]
        # get the x column:
        if args.x_column is not None:
            x_index = args.x_column
        else:
            x_index = 0
        # get the x data:
        x = data[:,x_index]
        if args.xlog: x = np.abs(x)
        # cycle over data columns:
        if args.y_columns is not None:
            data_indexes = args.y_columns
        else:
            data_indexes = range(0,n_col)    
        # remove the x axis:
        try:
            data_indexes.remove( x_index )
        except: pass
        # cycle over indexes:
        for ind in data_indexes:
            if args.ylog: 
                y = np.abs(data[:,ind])
            else:
                y = data[:,ind]
            
            plt.plot( x, y, label=str(ind) ) 
        
    plt.legend()
    
    # get log scales:
    if args.xlog: plt.xscale('log')
    if args.ylog: plt.yscale('log')
    
    # save or try to display:
    if args.outroot is not None:
        plt.savefig( args.outroot )
    else:
        plt.show()
    
    # close all:
    plt.close('all')
    
    # print some final feedback:
    if not args.quiet:
        if outroot is not None:
            print 'Done. Saved results in: ', outroot
        else:
            print 'Done.'

    # exit without error:
    exit(0)