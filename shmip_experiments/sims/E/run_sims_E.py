# -*- coding: utf-8 -*-
"""
SHMIP E simulations. 
"""

from dolfin import *
from constants import *
from sheet_model import *
from dolfin import MPI, mpi_comm_world
import time
import numpy as np 

#ns = range(1,6)
ns = [1]

MPI_rank = MPI.rank(mpi_comm_world())
# Input files 
input_files = ['../../inputs/E/input_E' + str(n) + '.hdf5' for n in ns]
# Result output directories
result_dirs = ['results_E' + str(n) for n in ns]

# Subdomain containing only a single outlet point at terminus
def outlet_boundary(x, on_boundary):
  cond1 = abs(x[0]) < 5.0
  cond2 = abs(x[1]) < 15.0
  return cond1 and cond2


for n in range(len(ns)):
  
  ### Setup the model  
  model_inputs = {}
  model_inputs['input_file'] = input_files[n]
  model_inputs['out_dir'] = result_dirs[n]
  model_inputs['constants'] = pcs
  # Point boundary condition at the outlet
  model_inputs['point_bc'] = outlet_boundary
  
  # Create the sheet model
  model = SheetModel(model_inputs)
  
  #f = interpolate(Constant(1.0), model.V_cg)
  #model.d_bcs[0].apply(f.vector())
  
    
  ### Run the simulation
  
  # Seconds per day
  spd = pcs['spd']
  # End time
  T = 1250.0 * spd
  # Time step
  dt = spd / 3.0
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
      model.write_pvds(['pfo', 'h', 'N'])
      
    if i % 8 == 0:
      model.checkpoint(['h', 'phi', 'N', 'q'])
    
    if MPI_rank == 0: 
      print
      
    i += 1
    
  end_time = time.time()
  np.savetxt(result_dirs[n] + '/Time_E' + str(ns[n]), np.array([start_time, end_time, end_time - start_time]))

  model.write_steady_file(result_dirs[n] + '/steady_E' + str(ns[n]))


  


 
 
