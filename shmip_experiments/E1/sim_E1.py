# -*- coding: utf-8 -*-
"""
SHMIP simulation E1. 
"""

from dolfin import *
from constants import *
from sheet_model import *
from dolfin import MPI, mpi_comm_world
from scale_functions import *

MPI_rank = MPI.rank(mpi_comm_world())
input_file = '../inputs/E1/inputs_E1.hdf5'
out_dir = 'results'


### Setup the model

model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['constants'] = pcs


# Create the sheet model
model = SheetModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 100.0 * spd
# Time step
dt = spd / 6.0
# Iteration count
i = 0

while model.t < T:  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)

  if i % 1 == 0:
    model.write_pvds(['pfo', 'h', 'q', 'phi'])
    
  if i % 1 == 0:
    model.checkpoint(['h', 'phi'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
  
model.write_steady_file(out_dir + '/steady')
