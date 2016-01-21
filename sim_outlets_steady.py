"""
Test simulation with two outlets with 0 pressure and the rest of the margin 0
flux. Runs to a steady state. 
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_outlets_steady'
# Input file
input_file = 'inputs/inputs_ref.hdf5'
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

# Use a smaller conductivity
pcs['k'] = 1e-2

model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['checkpoint_file'] = 'outlets_steady'
model_inputs['d_bcs'] = [bc]
model_inputs['constants'] = pcs


# Create the sheet model
model = SheetModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 75.0 * spd
# Time step
dt = 60.0 * 60.0 * 1.0
# Iteration count
i = 0

while model.t < T:
  if MPI_rank == 0: 
    current_time = model.t / spd
    #print "Current Time: " + str(current_time)
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds(['pfo', 'h'])
    
  if i % 1 == 0:
    model.checkpoint()
  
  if MPI_rank == 0: 
    print
    
  plot(model.pfo, interactive = True)
    
  i += 1