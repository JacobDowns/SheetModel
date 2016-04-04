# -*- coding: utf-8 -*-
"""
Steady state simulation. Trough bed. High melt. Spatially Varying k. 
"""

from dolfin import *
from constants import *
from channel_model import *
from dolfin import MPI, mpi_comm_world
from scale_functions import *

MPI_rank = MPI.rank(mpi_comm_world())


###  Create a spatially varying conductivity field

k_min = 5e-5
k_max = 5e-3

input_file = 'inputs_channel/inputs/inputs_trough_high.hdf5'

scale_functions = ScaleFunctions(input_file, k_min, k_max)
k = scale_functions.get_k(0.0)


### Setup the model

model_inputs = {}
pcs['k_c'] = 1e-1
model_inputs['k'] = k
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = 'channel_results/out_real_trough_steady_high/'
model_inputs['constants'] = pcs

# Create the sheet model
model = ChannelModel(model_inputs)

  
### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 550.0 * spd
# Time step
dt = 60.0 * 60.0 * 8.0
# Iteration count
i = 0

while model.t < T:  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds(['pfo', 'h', 'S'])
    
  if i % 1 == 0:
    model.checkpoint(['h', 'phi', 'S'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
  
model.write_steady_file('inputs_channel/steady/real_trough_steady_high')