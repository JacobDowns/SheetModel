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


# Process number
MPI_rank = MPI.rank(mpi_comm_world())


# Create a model mostly just to get the function space object
# Output directory
in_dir = 'sensitivity_test1/'
model_inputs = {}
model_inputs['input_file'] = 'inputs_sheet/inputs/inputs_high.hdf5'
model_inputs['out_dir'] = 'stuff'
model = SheetModel(model_inputs)
pfo = Function(model.V_cg)


# Load in the pfo field corresponding to different k, ub combinations
ks = linspace(1e-5, 5e-3, 40)
ubs = linspace(0, 150, 40)

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

kg = kg.reshape((40, 40))
ug = ug.reshape((40, 40))
pfos1 = pfos1.reshape((40, 40))
pfos2 = pfos2.reshape((40, 40))
pfos3 = pfos3.reshape((40, 40))


# Make contour plots of pressures at the three test points

pfo_min = min([pfos1.min(), pfos2.min(), pfos3.min()])
pfo_max = min([pfos1.max(), pfos2.max(), pfos3.max()])

subplot(1,3,1)

levels = arange(0, 1, 0.005)
CS1 = contourf(kg, ug, pfos1, levels = levels, cmap=plt.cm.jet, vmin = pfo_min, vmax = pfo_max)
CS2 = contour(kg, ug, pfos1, levels = [0.96, 0.94, 0.92, 0.9, 0.88, 0.86, 0.84, 0.82, 0.8, 0.78, 0.76, 0.74, 0.72] , colors='k')
plot(kg, ug, 'ko', ms=1)

subplot(1,3,2)
CS3 = contourf(kg, ug, pfos2, 100, cmap=plt.cm.jet, vmin = pfo_min, vmax = pfo_max)
CS4 = contour(kg, ug, pfos2, 20, colors='k')
plot(kg, ug, 'ko', ms=1)
colorbar(CS3)

subplot(1,3,3)
CS5 = contourf(kg, ug, pfos3, 100, cmap=plt.cm.jet, vmin = pfo_min, vmax = pfo_max)
CS6 = contour(kg, ug, pfos3, 20, colors='k')
plot(kg, ug, 'ko', ms=1)

#clabel(CS, inline=1, fontsize=16, fmt='%1.2f')
#clabel(CS, inline=1, fontsize=10)
show()