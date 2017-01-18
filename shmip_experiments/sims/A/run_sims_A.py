# -*- coding: utf-8 -*-
"""
SHMIP A simulations. 
"""

from dolfin import *
from constants import *
from sheet_model import *
from dolfin import MPI, mpi_comm_world
import time
import numpy as np 

ns = [1]

MPI_rank = MPI.rank(mpi_comm_world())
input_files = ['../../inputs/A/input_A' + str(n) + '.hdf5' for n in ns]
result_dirs = ['resultsc_A' + str(n) + '_c' for n in ns]

for n in range(len(ns)):
  
  ### Setup the model
  
  model_inputs = {}
  model_inputs['input_file'] = input_files[n]
  model_inputs['out_dir'] = result_dirs[n]
  model_inputs['constants'] = pcs

  # Create the sheet model
  model = SheetModel(model_inputs)
  
  #k = interpolate(Constant(9e-3), model.V_cg)
  model.set_k(k)
  

  ### Run the simulation
  
  # Seconds per day
  spd = pcs['spd']
  # End time
  T = 1500.0 * spd
  # Time step
  dt = spd / 4.0
  # Iteration count
  i = 0
  
  # Time the run  
  start_time = time.time()
  
  while model.t < T:  
    if MPI_rank == 0: 
      current_time = model.t / spd
      print 'Current time: ' + str(current_time)
    
    model.step(dt)
    
    if i % 8 == 0:
      model.write_pvds(['pfo', 'h'])
      
    if i % 4 == 0:
      model.checkpoint(['h', 'phi', 'N'])
    
    if MPI_rank == 0: 
      print
      
    i += 1
  

  model.write_steady_file(result_dirs[n] + '/steady_A' + str(ns[n]))

  end_time = time.time()
  np.savetxt(result_dirs[n] + '/Time_A' + str(ns[n]), np.array([start_time, end_time, end_time - start_time]))
  


 
 
