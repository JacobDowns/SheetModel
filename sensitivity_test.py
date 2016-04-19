"""
Solve for pressure with a bunch of constant k, ub combinations with sheet 
height fixed. 
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *
from numpy import *
import itertools


# Process number
MPI_rank = MPI.rank(mpi_comm_world())

# Output directory
out_dir = 'sensitivity_test'

model_inputs = {}
model_inputs['input_file'] = 'inputs_sheet_/steady/ref_steady.hdf5'
model_inputs['out_dir'] = 'sensitivity_test/'

model_inputs['opt_params']= {'tol' : 1e-3, 'scale' : 50}

ks = linspace(5e-5, 5e-3, 25)
ubs = linspace(0, 200, 25) / pcs['spy']

model = SheetModel(model_inputs)
h = Function(model.V_cg)
h.assign(model.h)
model.set_m(Function(model.V_cg))

xs = project(Expression("x[0]"), model.V_cg)
plot(xs, interactive = True)
quit()

i = 0
for c in itertools.product(ks, ubs):
  k = c[0]
  ub = c[1]
  
  if MPI_rank == 0:
    print (i, len(ks) * len(ubs))
    print c
    print ("k", k, "ub", ub)

  model.set_h(h)
  model.set_k(interpolate(Constant(k), model.V_cg))
  model.set_u_b(interpolate(Constant(ub), model.V_cg)) 
  
  model.step(1.0)
  
  File(out_dir + "/pfo_" + str(i) + ".xml") << model.pfo
  i += 1
  

