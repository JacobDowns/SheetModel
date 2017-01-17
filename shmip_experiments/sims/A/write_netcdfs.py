# -*- coding: utf-8 -*-
"""
Write results as netcdf files
"""

from steady_view import *

ns = [1,2,3]

for n in ns:
  # Load steady state results file
  view = SteadyView('results_A' + str(n) + '/steady_A' + str(n) + '.hdf5')
  # Write results as netcdf file
  view.write_netcdf('results_netcdf/A' + str(n), 'Jacob Z Downs, A' + str(n))