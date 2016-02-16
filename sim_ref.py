"""
Test simulation to make sure nothing is broken.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


### Model inputs

# Process number
MPI_rank = MPI.rank(mpi_comm_world())

# File with model inputs for the reference experiment

model_inputs = {}
model_inputs['input_file'] = 'out_ref/out.hdf5'
#model_inputs['input_file'] = 'inputs/inputs_ref.hdf5'
model_inputs['out_dir'] = 'out_ref/'

# Create the sheet model
model = SheetModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 75.0 * spd
# Time step
dt = spd
# Iteration count
i = 0

while model.t < T:  
  # Update the melt
  model.set_m(project(Constant((model.t / spd) * 1e-10), model.V_cg))
  
  if MPI_rank == 0: 
    current_time = model.t / spd
    #print "Current Time: " + str(current_time)
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds(['pfo', 'h'])
    
  if i % 1 == 0:
    model.checkpoint(['m'])
  
  if MPI_rank == 0: 
    print
    
  i += 1