# -*- coding: utf-8 -*-
"""
Run the channel model to steady state with a variety of k, k_c combinations.
"""

from dolfin import *
from constants import *
from channel_model import *
from dolfin import MPI, mpi_comm_world
import itertools
import numpy as np

MPI_rank = MPI.rank(mpi_comm_world())

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

model_inputs = {}
model_inputs['input_file'] = 'inputs_channel/inputs_channel_ref.hdf5'
model_inputs['newton_params'] = prm

k_min = 5e-4
k_max = 5e-2
k_c_min = 5e-2
k_c_max = 5e-1
ks = np.linspace(k_min, k_max, 5)
k_cs = np.linspace(k_c_min, k_c_max, 5)


# Run a model to steady state
def run():  
  # Seconds per day
  spd = pcs['spd']
  # End time
  T = 500.0 * spd
  # Time step
  dt = spd / 2.0
  # Iteration count
  i = 0
  
  while model.t < T:  
    
    if MPI_rank == 0: 
      current_time = model.t / spd
      #print "Current Time: " + str(current_time)
      print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
    
    model.step(dt)
    
    if i % 1 == 0:
      model.write_pvds(['pfo', 'h', 'S'])
      
    if i % 1 == 0:
      model.checkpoint(['h', 'phi', 'S'])
    
    if MPI_rank == 0: 
      print
      
    i += 1


# Create a bunch of models with different conductivities and run them
for c in itertools.product(ks, k_cs):
  k = c[0]
  k_c = c[1]
  pcs['k'] = k
  pcs['k_c'] = k_c
  model_inputs['constants'] = pcs
  model_inputs['out_dir'] = 'out_channel_k_' + str(k) + "_k_c_" + str(k_c)
  model_inputs['checkpoint_file'] = 'out_channel_k_' + str(k) + "_k_c_" + str(k_c)

  model = ChannelModel(model_inputs)
  run()

