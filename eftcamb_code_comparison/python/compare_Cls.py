import matplotlib
matplotlib.use('Agg')

import CAMB_plots_lib.CAMB_comp_plots   as Cplt
import CAMB_plots_lib.plot_colors       as pltcol

import matplotlib.pyplot as plt
import numpy as np
import math
import sys

# parse command line arguments
root1    = sys.argv[1]
root2    = sys.argv[2]
outroot  = sys.argv[3]

print "Comparing: ", root1, " and ", root2

if len(sys.argv) == 6:
    label1   = sys.argv[4]
    label2   = sys.argv[5]

    result = Cplt.CAMB_results_compare_plot(root1, root2, outroot,
                                        lensing=True, transfer=True, tensor=True,
                                        name1=label1, name2=label2)
else:
    result = Cplt.CAMB_results_compare_plot(root1, root2, outroot,
                                        lensing=True, transfer=True, tensor=True)

# change the colors
result.color1      = pltcol.plot_colors.color_scheme_1D.default[0]
result.color2      = pltcol.plot_colors.color_scheme_1D.default[1]
result.color_compa = pltcol.plot_colors.color_scheme_1D.default[3]

# make the plots:
try:
    result.plot_compare_scalCls()
except Exception,e:
    print 'Problem plotting Scalar Cls:'
    print e
try:
    result.plot_compare_lensedCls()
except Exception,e:
    print 'Problem plotting Lensed Cls:'
    print e
try:
    result.plot_compare_tensCls()
except Exception,e:
    print 'Problem plotting Tensor Cls:'
    print e
try:
    result.plot_compare_totalCls()
except Exception,e:
    print 'Problem plotting Total Cls:'
    print e
try:
    result.plot_compare_scalCovCls()
except Exception,e:
    print 'Problem plotting Scalar Cls Cov:'
    print e
try:
    result.plot_compare_Transfer()
except Exception,e:
    print 'Problem plotting Transfer:'
    print e

# exit if the file is directly called:
if __name__ == "__main__":
    exit()
