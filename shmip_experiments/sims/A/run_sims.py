# -*- coding: utf-8 -*-
"""
SHMIP A simulations. 
"""

from dolfin import *
from constants import *
from sheet_model import *
from dolfin import MPI, mpi_comm_world

ns = range(1,7)

MPI_rank = MPI.rank(mpi_comm_world())
input_files = '../inputs/A1/inputs_A1.hdf5'

input_files = ['../inputs/A' + str(n) + '.hdf5' for n in ns]
result_dirs = ['results_A' + str(n) for n in ns]

print input_files
print result_dirs
quit()


### Setup the model

model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['constants'] = pcs

# Create the sheet model
model = ChannelModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 5000.0 * spd
# Time step
dt = spd / 4.0
# Iteration count
i = 0

while model.t < T:  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print 'Current time: ' + str(current_time)
  
  model.step(dt)
  
  if i % 8 == 0:
    model.write_pvds(['pfo', 'h'])
    
  if i % 4 == 0:
    model.checkpoint(['h', 'phi', 'S', 'k', 'N'])
  
  if MPI_rank == 0: 
    print
    
  i += 1

model.write_steady_file(out_dir + '/steady')