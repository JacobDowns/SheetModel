"""
Runs Issunguata to a steady state. Zero bcs.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_is_steady'
# Input file
input_file = 'inputs/inputs_is.hdf5'


### Initialize model

pcs['k'] = 7e-3
model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['opt_params'] = {'tol' : 1e-2, 'scale' : 30}

# Create the sheet model
model = SheetModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 90.0 * spd
# Time step
dt = 60.0 * 60.0 * 1.0
# Iteration count
i = 0

while model.t < T:
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds(['pfo', 'h'])
    
  if i % 1 == 0:
    model.checkpoint()
  
  if MPI_rank == 0: 
    print
    
  i += 1