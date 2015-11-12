"""
Test simulation to make sure nothing is broken.
"""

import sys
from dolfin import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '/home/jake/schoof/')
from sheet_model import *
from constants import *

# Model input directory
in_dir = "inputs_high_melt/"
# Output directory
out_dir = "out_test/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())

# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)


### Model inputs

# Initial sheet height
h_init = Function(V_cg)
h_init.interpolate(Constant(0.05))

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['out_dir'] = out_dir

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

