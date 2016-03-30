# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 16:02:10 2016

@author: jake
"""

from dolfin import *
from constants import *
from channel_model import *
from dolfin import MPI, mpi_comm_world

MPI_rank = MPI.rank(mpi_comm_world())

# Output directory
out_dir = 'out_channel_outlets'
# Input file
input_file = 'inputs_channel/inputs_channel_ref.hdf5'
# Load the input file
inputs = HDF5File(mpi_comm_world(), input_file, 'r')


### Create outlet boundary conditions
mesh = Mesh()
inputs.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1)

def outlet_boundary(x, on_boundary):
  cond1 = (abs(x[0] - 150.0) < 150.0) and (abs(x[1] - 6666.0) < 150.0)
  cond2 = (abs(x[0] - 150.0) < 150.0) and (abs(x[1] - 13332.0) < 150.0)
  return cond1 or cond2

bc = DirichletBC(V_cg, 0.0, outlet_boundary, "pointwise")


### Initialize model
prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['constants'] = pcs
model_inputs['newton_params'] = prm
model_inputs['d_bcs'] = [bc]

# Create the sheet model
model = ChannelModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 300.0 * spd
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