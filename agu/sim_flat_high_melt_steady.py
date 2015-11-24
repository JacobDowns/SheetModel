"""
Steadu state fpr the high melt scenario. Distributed melt.
"""

import sys
from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '../')
from sheet_model import *
from constants import *

# Model input directory
in_dir = "inputs_flat_bed/"
# Output directory
out_dir = "out_flat_bed_high_melt_steady/"
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

# Distributed melt
m = Function(V_cg)
File(in_dir + "m.xml") >> m

# Minimum conductivity
k_min = 5e-5
# Scaling parameter
a = 31220.5955507
# Hydraulic conductivity
k = project(a * m + k_min, V_cg)

# Newton solver params
prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-5
prm['newton_solver']['absolute_tolerance'] = 1e-5
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm
model_inputs['m'] = m
model_inputs['k'] = k
model_inputs['constants'] = pcs

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 100.0 * spd
# Time step
dt = 60.0 * 60.0 * 12.0
# Iteration count
i = 0

while model.t < T:  
  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 2 == 0:
    model.write_pvds(['h', 'u_b', 'm', 'pfo', 'k'])
    
  if i % 2 == 0:
    model.write_xmls(['h', 'pfo', 'phi'])
  
  if MPI_rank == 0: 
    print
    
  i += 1

