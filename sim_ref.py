"""
Test simulation to make sure nothing is broken.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *



### Model inputs

# File with model inputs for the reference experiment
input_file = "inputs/inputs_ref.hdf5"

# Initial sheet height
h_init = Function(V_cg)
h_init.interpolate(Constant(0.05))

model_inputs = {}
model_inputs['input_file'] = input_file

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 75.0 * spd
# Time step
dt = 60.0 * 60.0 * 8.0
# Irataion count
i = 0

while model.t < T:
  if MPI_rank == 0: 
    current_time = model.t / spd
    #print "Current Time: " + str(current_time)
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds()
    
  if i % 1 == 0:
    model.write_xmls()
  
  if MPI_rank == 0: 
    print
    
  i += 1


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 75.0 * spd
# Time step
dt = 60.0 * 60.0 * 8.0
# Irataion count
i = 0

while model.t < T:
  if MPI_rank == 0: 
    current_time = model.t / spd
    print "Current Time: " + str(current_time)
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds()
    
  if i % 1 == 0:
    model.write_xmls()
  
  if MPI_rank == 0: 
    print
    
  i += 1