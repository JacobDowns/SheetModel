"""
Steady state for the zero bc simulation with spatially varying conductivity.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_is_steady_realistic'
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
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['constants'] = pcs
model_inputs['opt_params'] = {'tol' : 5e-3, 'scale' : 30}
model = SheetModel(model_inputs)


### Make the conductivity spatially varying depending on melt
m = Function(model.V_cg)
m.assign(model.m)
# Minimum conductivity
k_min = 9e-5
# Scaling parameter that sets the maximum possible conductivity
a = model.pcs['k'] / m.vector().max()

model.set_k(Constant(a)*m + Constant(k_min))
plot(k, interactive = True)
quit()


### Run the simulation

# Seconds per day
spd = pcs['spd']
# Seconds per month
spm = pcs['spm']
# End time
T = 7.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Iteration count
i = 0


while model.t < T:  
  # Update the conductivity
  model.set_k(project(Constant(k_scale(model.t)) * m + Constant(k_min), model.V_cg))

  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 3 == 0:
    model.write_pvds(['h', 'pfo'])
    
  if i % 3 == 0:
    model.checkpoint(['pfo', 'h', 'k', 'm', 'u_b'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
