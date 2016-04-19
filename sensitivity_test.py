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
out_dir = 'paper_results/sensitivity_test'

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 35

model_inputs = {}
model_inputs['input_file'] = 'inputs_sheet/steady/ref_steady.hdf5'
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm

ks = linspace(5e-5, 5e-3, 30)
ubs = linspace(0, 100, 30) / pcs['spy']

model = SheetModel(model_inputs)
h = Function(model.V_cg)
h.assign(model.h)
model.set_m(Function(model.V_cg))

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
  
  model.checkpoint(['k', 'u_b', 'pfo'])
  
  File(out_dir + "/pfo_" + str(i) + ".xml") << model.pfo
  i += 1
  

