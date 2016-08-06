import CAMB_plots_lib.CAMB_single_plots   as Cplt
import CAMB_plots_lib.plot_colors         as pltcol

import matplotlib.pyplot as plt
import numpy as np
import math
import sys

# parse command line arguments
root    = sys.argv[1]
outroot = sys.argv[2]

if len(sys.argv) == 3:
    label  = sys.argv[3]
    result = Cplt.CAMB_results_plot(root, outroot, lensing=True, transfer=True, tensor=True, name=label)
else:
    result = Cplt.CAMB_results_plot(root, outroot, lensing=True, transfer=True, tensor=True)

# change the colors
result.color = pltcol.plot_colors.color_scheme_1D.default[0]

# make the plots
result.plot_scalCls()
result.plot_lensedCls()
result.plot_tensCls()
result.plot_totalCls()
result.plot_scalCovCls()
result.plot_Transfer()