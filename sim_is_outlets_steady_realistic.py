"""
Steady state for the outlet simulation with spatially varying conductivity.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_is_outlets_steady_realistic'
# Input file
input_file = 'inputs/inputs_is_outlets_steady.hdf5'
# Load the input file
inputs = HDF5File(mpi_comm_world(), input_file, 'r')


### Create outlet boundary conditions
mesh = Mesh()
inputs.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1)

# Load potential at 0 pressure
B = Function(V_cg)
inputs.read(B, "B")

def outlet_boundary(x, on_boundary):
  # These two outlet points are based on Google Earth -- places where it looks
  # like water is flowing out
  out1_x = -491254.66
  out1_y = -2461159.0
  out2_x = -491839.3
  out2_y = -2472998.0
  
  cond1 = (abs(x[0] - out1_x) < 150.0) and (abs(x[1] - out1_y) < 150.0)
  cond2 = (abs(x[0] - out2_x) < 150.0) and (abs(x[1] - out2_y) < 150.0)
  return cond1 or cond2

bc = DirichletBC(V_cg, pcs['rho_w'] * pcs['g'] * B, outlet_boundary, "pointwise")


### Initialize model

# Use a smaller conductivity
pcs['k'] = 7e-3
model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['d_bcs'] = [bc]
model_inputs['constants'] = pcs
model_inputs['opt_params'] = {'tol' : 1e-2, 'scale' : 30}
model = SheetModel(model_inputs)


### Make the conductivity spatially varying depending on melt
m = Function(model.V_cg)
m.assign(model.m)

# Maximum conductivity
k_max = model.pcs['k']
# Minimum conductivity
k_min = 7e-5
# Scaling parameter that sets the maximum possible conductivity
a = (k_max - k_min) / m.vector().max()
# Set the conductivity
model.set_k( project(Constant(a)*m + Constant(k_min), V_cg))


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 90 * spd
# Time step
dt = spd
# Iteration count
i = 0


while model.t < T:  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  plot(model.pfo)
  quit()
  
  if i % 1 == 0:
    model.write_pvds(['h', 'pfo'])
    
  if i % 1 == 0:
    model.checkpoint(['h'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
