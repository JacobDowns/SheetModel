# -*- coding: utf-8 -*-
"""
Steady state sim. Spatially dependent conductivity. Trough. 
"""

from dolfin import *
from constants import *
from channel_model import *
from dolfin import MPI, mpi_comm_world

MPI_rank = MPI.rank(mpi_comm_world())

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

model_inputs = {}
pcs['k'] = 5e-3
pcs['k_c'] = 1e-1
model_inputs['input_file'] = 'inputs_channel/inputs_channel_ref_trough.hdf5'
model_inputs['out_dir'] = 'channel_realistic/out_channel_ref_trough_realistic/'
model_inputs['constants'] = pcs
model_inputs['newton_params'] = prm

# Create the sheet model
model = ChannelModel(model_inputs)


### Compute the conductivity

m = Function(model.V_cg)
m.assign(model.m)

# Maximum conductivity
k_max = model.pcs['k']
# Minimum conductivity
k_min = 5e-5
# Scaling parameter that sets the maximum possible conductivity
a = (k_max - k_min) / m.vector().max()
# Set the conductivity
model.set_k( project(Constant(a)*m + Constant(k_min), model.V_cg))
  

### Run the simulation

# Seconds per day
spd = pcs['spd']
# Seconds per month
spm = pcs['spm']
# End time
T = 550.0 * spd
# Time step
dt = 60.0 * 60.0 * 8.0
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