# -*- coding: utf-8 -*-
"""
Gammut of low melt lag time tests. Week shutoff.
"""

from dolfin import *
from constants import *
from channel_model import *
from dolfin import MPI, mpi_comm_world
from scale_functions import *

MPI_rank = MPI.rank(mpi_comm_world())


### Load the input file 

input_file = 'inputs_channel/steady/real_steady_low.hdf5'
k_min = 5e-5
k_max = 5e-3
scale_functions = ScaleFunctions(input_file, k_min, k_max, u_b_max = 80.0, shutoff_time = 7.0 * pcs['spd'])


### Setup the model

# Seconds per day
spd = pcs['spd']
model_inputs = {}
pcs['k_c'] = 1e-1
model_inputs['constants'] = pcs
model_inputs['input_file'] = input_file

lag_times = [0.0, spd, 7.0 * spd, 30.0 * spd]
names = ['no', 'day', 'week', 'month']


### Run the simulation

def run_sim():

  # Seconds per month
  spm = pcs['spm']
  # End time
  T = 8.0 * spm
  # Time step
  dt = 60.0 * 60.0 * 8.0
  # Iteration count
  i = 0
  
  while model.t < T:  
    # Update the melt
    model.set_m(scale_functions.get_m(model.t))
    # Update the conductivity
    model.set_k(scale_functions.get_k(model.t))
    # Update the sliding speed
    model.set_u_b(scale_functions.get_u_b(model.t))  
    
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
    
    
### Run simulations with different lag times 

i = 0
for b in lag_times:
  scale_functions.b = b
  name = names[i]
  model_inputs['out_dir'] = 'channel_results/real_lag_low1/' +  name + '_lag'
  model_inputs['checkpoint_file'] = 'checkpoint/' +  name + '_lag'
  model = ChannelModel(model_inputs)
  run_sim()
  i += 1
