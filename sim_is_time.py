"""
IS winter pressure simulation with fixed conductivity and sliding.
"""

from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *

# Model input directory
in_dir = "inputs/"
# Output directory
out_dir = "out_is_winter/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())


### Model inputs

model_inputs = {}
# Load the steady state IS input file
model_inputs['input_file'] = 'inputs/inputs_is_steady.hdf5'
model_inputs['out_dir'] = 'out_is_winter'

# Create the sheet model
model = SheetModel(model_inputs)


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