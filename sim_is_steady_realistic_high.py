"""
Steady state for the zero bc simulation with spatially varying conductivity. 
1.5x the normal melt.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_is_steady_realistic_high'
# Input file
input_file = 'inputs/inputs_is_steady.hdf5'
# Load the input file
inputs = HDF5File(mpi_comm_world(), input_file, 'r')


### Create outlet boundary conditions
mesh = Mesh()
inputs.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1)


### Initialize model

model_inputs = {}
pcs['k'] = 1e-2
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['constants'] = pcs
model_inputs['opt_params'] = {'tol' : 1e-2, 'scale' : 30}
model = SheetModel(model_inputs)


### Make the conductivity spatially varying depending on melt
model.set_m(project(model.m * 1.5, V_cg))
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
  
  if i % 1 == 0:
    model.write_pvds(['h', 'pfo'])
    
  if i % 1 == 0:
    model.checkpoint(['h'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
