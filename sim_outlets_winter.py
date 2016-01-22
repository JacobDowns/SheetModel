"""
Starts from the steady state generated by sim_outlets_steady and reduces melt
while keeping sliding and conductivity constant.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_outlets_winter'
# Input file
input_file = 'inputs/inputs_outlets_steady.hdf5'
# Load the input file
inputs = HDF5File(mpi_comm_world(), input_file, 'r')


### Create outlet boundary conditions
mesh = Mesh()
inputs.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1)

def outlet_boundary(x, on_boundary):
  cond1 = (abs(x[0] - 150.0) < 150.0) and (abs(x[1] - 6666.0) < 150.0)
  cond2 = (abs(x[0] - 150.0) < 150.0) and (abs(x[1] - 13332.0) < 150.0)
  return cond1 or cond2

bc = DirichletBC(V_cg, 0.0, outlet_boundary, "pointwise")


### Melt scale function 

# Copy the summer high melt function
m = Function(model.V_cg)
m.assign(model.m)

# Seconds per day
spd = pcs['spd']
# How long it takes for the melt to shut off completely
shutoff_length = 30.0 * spd
 
# Function that reduces melt to 0 over time
def m_scale(t):
  if t < 0.0:
    return 1.0
  if t <= shutoff_length:
    return cos((pi / (2.0 * shutoff_length)) * t)
  return 0.0


### Initialize model

# Use a smaller conductivity
pcs['k'] = 1e-2

model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['checkpoint_file'] = 'outlets_steady'
model_inputs['d_bcs'] = [bc]
model_inputs['constants'] = pcs

# Create the sheet model
model = SheetModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# Seconds per month
spm = pcs['spm']
# End time
T = 10.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Iteration count
i = 0

while model.t < T:  
  
  # Update the melt
  model.set_m(project(Constant(m_scale(model.t)) * m, model.V_cg))
  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 3 == 0:
    model.write_pvds(['h', 'pfo'])
    
  if i % 3 == 0:
    model.checkpoint(['pfo'])
  
  if MPI_rank == 0: 
    print
    
  i += 1