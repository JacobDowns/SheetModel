# -*- coding: utf-8 -*-
"""
Winter sim. Time and spatially varying k. Flat bed. 
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
model_inputs['input_file'] = 'inputs_channel/channel_ref_steady_realistic.hdf5'
model_inputs['out_dir'] = 'channel_realistic/out_channel_ref_realistic_winter_week_lag/'
model_inputs['constants'] = pcs
model_inputs['newton_params'] = prm

# Create the sheet model
model = ChannelModel(model_inputs)


### Set up scaling functions for time dependent m, k, and u_b

m = Function(model.V_cg)
m.assign(model.m)
u_b = Function(model.V_cg)
u_b.assign(model.u_b)

# Seconds per day
spd = pcs['spd']
# Seconds per year
spy = pcs['spy']
# How long it takes for the melt to shut off completely
shutoff_length = 30.0 * spd
 
# Function that reduces melt to 0 over time
def m_scale(t):
  if t < 0.0:
    return 1.0
  if t <= shutoff_length:
    return cos((pi / (2.0 * shutoff_length)) * t)
  return 0.0

# Maximum conductivity
k_max = model.pcs['k']
# Minimum conductivity
k_min = 5e-5
# Scaling parameter that sets the maximum possible conductivity
a = (k_max - k_min) / m.vector().max()
# Lag 
b = 7.0 * spd
# Function that scales k proportionally to m
def k_scale(t):
  return a * m_scale(t - spd)
  
# Create a function that scales u_b down to a maximum of 80 (m/a) in winter
c = (80.0 / pcs['spy']) / u_b.vector().max()
def u_b_scale(t):
  return -m_scale(t) * (c - 1.0) + c
  

### Run the simulation

# Seconds per day
spd = pcs['spd']
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
  model.set_m(project(Constant(m_scale(model.t)) * m, model.V_cg))
  # Update the conductivity
  model.set_k(project(Constant(k_scale(model.t)) * m + Constant(k_min), model.V_cg))
  # Update the sliding speed
  model.set_u_b(project(Constant(u_b_scale(model.t)) * u_b, model.V_cg))  
  
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