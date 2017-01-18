# -*- coding: utf-8 -*-
"""
SHMIP D simulations. 
"""

from dolfin import *
from constants import *
from sheet_model import *
from dolfin import MPI, mpi_comm_world
import time
import numpy as np 

ns = [1]

MPI_rank = MPI.rank(mpi_comm_world())
# Input file is steady state from A1
input_file = '../../inputs/D/steady_A1.hdf5'
#input_file = 'results_D1/year7.hdf5'
# Result output directories
result_dirs = ['results_D' + str(n) for n in ns]

## Params

# Delta temps for each run
DTs = [-4, -2, 0, 2, 4]
spy = pcs['spy']
# Lapse rate (K / m)
lr = -0.0075 
# Degre day factor (m/K/s)
DDF = 0.01 / 86400.0 

for n in range(len(ns)):
  
  ### Setup the model
  
  model_inputs = {}
  model_inputs['input_file'] = input_file
  model_inputs['out_dir'] = result_dirs[n]
  model_inputs['constants'] = pcs
  model_inputs['checkpoint_file'] = 'out1'

  # Create the sheet model
  model = SheetModel(model_inputs)

  # Temperature at 0m  
  def temp(t):
    return -16.0 * np.cos(2.0*pi/spy*t) - 5.0 + DTs[n]

  # Melt
  m = Function(model.V_cg)
  # Array of zeros, same length as local array in m
  zs = np.zeros(len(m.vector().array()))
  # Calculate melt
  def update_m(t):    
    ms = np.maximum(zs, (model.H.vector().array()*lr + temp(t))*DDF) + 7.93e-11
    m.vector().set_local(ms)
    m.vector().apply("insert")

    
  ### Run the simulation
  
  # Seconds per day
  spd = pcs['spd']
  # Seconds per year
  spy = pcs['spy']
  # Time step
  dt = spd / 20.0
  # Iteration count
  i = 0
  # Time the run  
  start_time = time.time()
  # Number of years
  year = 1
  

  # Put output for each year in a separate folder
  while year <= 15:
    
    pvd_dir = result_dirs[n] + '/year' + str(year) + '/'
    out_pfo = File(pvd_dir + 'pfo.pvd')
    out_h = File(pvd_dir + 'h.pvd')
    out_N = File(pvd_dir + 'N.pvd')
    
    # Run the model for a year
    while model.t < year * spy:  
      if MPI_rank == 0: 
        current_time = model.t / spd
        print 'Current time: ' + str(current_time)
        
      # Update melt
      update_m(model.t)
      model.set_m(m)
      
      model.step(dt)
      
      if i % 20 == 0:
        out_pfo << model.pfo
        out_h << model.h
        out_N << model.N
        
      if i % 20 == 0:
        model.checkpoint(['h', 'phi', 'N', 'm', 'q'])
      
      if MPI_rank == 0: 
        print
        
      i += 1
    
    end_time = time.time()
    np.savetxt(result_dirs[n] + '/Time_year' + str(year), np.array([start_time, end_time, end_time - start_time]))
    
    model.write_steady_file(result_dirs[n] + '/year' + str(year))
    year += 1


  


 
 
