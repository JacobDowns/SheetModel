"""
Reference simulation over winter.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *
from pylab import *
import itertools


# Process number
MPI_rank = MPI.rank(mpi_comm_world())

model_inputs = {}
model_inputs['input_file'] = 'inputs/steady_ref/steady_ref.hdf5'
out_dir = 'sensitivity_test'
model_inputs['out_dir'] = 'sensitivity_test/'

ks = linspace(5e-6, 5e-3, 25)
ubs = linspace(0, 200, 25)

d = {}

i = 0
for c in itertools.product(ks, ubs):
  print (i, len(ks) * len(ubs))
  k = c[0]
  ub = ubs[0]
  
  print c
  print ("k", k, "ub", ub)
  
  model_inputs['checkpoint_file'] = 'k_' + str(k) + "_ub_" + str(ub)
  model = SheetModel(model_inputs)
  
  model.set_m(Function(model.V_cg))
  model.set_k(interpolate(Constant(k), model.V_cg))
  model.set_u_b(interpolate(Constant(ub), model.V_cg)) 
  
  model.step(1.0)
  
  File(out_dir + "/pfo_" + str(i) + ".xml") << model.pfo
  i += 1
