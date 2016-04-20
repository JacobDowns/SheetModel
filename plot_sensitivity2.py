# -*- coding: utf-8 -*-
"""
Plot the results of the sensitivity test.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *
from pylab import *
import itertools
import pickle
from matplotlib.colors import LinearSegmentedColormap


# Process number
MPI_rank = MPI.rank(mpi_comm_world())

matplotlib.rcParams.update({'font.size': 20})


# Create a model mostly just to get the function space object
# Output directory
in_dir = 'sensitivity_test/'
model_inputs = {}
model_inputs['input_file'] = 'inputs_sheet/inputs/inputs_high.hdf5'
model_inputs['out_dir'] = 'stuff'
model = SheetModel(model_inputs)
pfo = Function(model.V_cg)


# Load in the pfo field corresponding to different k, ub combinations
ks = linspace(1e-5, 5e-3, 30)
ubs = linspace(0, 100, 30)

i = 0

kg = []
ug = []
pfos1 = []
pfos2 = []
pfos3 = []

for c in itertools.product(ks, ubs):
  print i
  k = c[0]
  ub = c[1]
  
  kg.append(k)
  ug.append(ub)
  
  # Load in pressure
  File(in_dir + "pfo_" + str(i) + ".xml") >> pfo
  
  #dolfin.plot(pfo, interactive = True)
  pfos1.append(pfo([10e3, 10e3]))  
  pfos2.append(pfo([20e3, 10e3]))
  pfos3.append(pfo([50e3, 10e3]))
  
  i += 1
 
kg = array(kg)
ug = array(ug)
pfos1 = array(pfos1)
pfos2 = array(pfos2)
pfos3 = array(pfos3)

kg = kg.reshape((30, 30))
ug = ug.reshape((30, 30))
pfos1 = pfos1.reshape((30, 30))
pfos2 = pfos2.reshape((30, 30))
pfos3 = pfos3.reshape((30, 30))


# Make contour plots of pressures at the three test points

pfo_min = min([pfos1.min(), pfos2.min(), pfos3.min()])
pfo_max = min([pfos1.max(), pfos2.max(), pfos3.max()]) + 1e-2



cdict = pickle.load(open("sens.p", "rb"))
desaturated = LinearSegmentedColormap('desaturated', cdict)


levels = arange(pfo_min - (pfo_min % 0.005), pfo_max, 0.005)
cf_levels1 = arange(pfos1.min() - (pfos1.min() % 0.02), pfos1.max(), 0.02)
cf_levels2 = arange(pfos2.min() - (pfos2.min() % 0.02), pfos2.max(), 0.02)
cf_levels3 = arange(pfos2.min() - (pfos2.min() % 0.02), pfos3.max(), 0.02)


subplot(1,3,1)
title('1')
CS1 = contourf(kg, ug, pfos1, levels = levels, cmap=desaturated, vmin = pfo_min, vmax = pfo_max)
CS2 = contour(kg, ug, pfos1, levels = cf_levels1 , colors='k')
clabel(CS2, levels = cf_levels1, inline=1, fmt='%0.2f', fontsize=18)
ylabel(r'Sliding Speed (ma$^{-1}$)')
plot(kg, ug, 'ko', ms=1)

subplot(1,3,2)
title('2')
CS3 = contourf(kg, ug, pfos2, levels = levels, cmap=desaturated, vmin = pfo_min, vmax = pfo_max)
CS4 = contour(kg, ug, pfos2, levels = cf_levels2, colors='k')
clabel(CS4, levels = cf_levels2, inline=1, fmt='%1.2f', fontsize=18)
xlabel(r'Hydraulic Conductivity (m$^{7/4}$kg$^{-1/2}$)')
plot(kg, ug, 'ko', ms=1)

subplot(1,3,3)
title('3')
CS5 = contourf(kg, ug, pfos3, levels = levels, cmap=desaturated, vmin = pfo_min, vmax = pfo_max)
CS6 = contour(kg, ug, pfos3, cf_levels3, colors='k')
clabel(CS6, levels = cf_levels3, inline=1, fmt='%0.2f', fontsize=18)
plot(kg, ug, 'ko', ms=1)


colorbar(CS3, ticks=arange(0.74, 1.0, 0.02), label = "PFO")
#clabel(CS, inline=1, fontsize=16, fmt='%1.2f')

show()