"""
IS steady state simulation.
"""

import sys
from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '../')
from sheet_model import *
from constants import *

# Model input directory
in_dir = "inputs_is/"
# Output directory
out_dir = "out_is_steady/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)


### Model inputs

# Initial sheet height
h_init = Function(V_cg)
h_init.interpolate(Constant(0.05))

# Sliding speed
u_b = Function(V_cg)
File(in_dir + "u_b.xml") >> u_b

m = Function(V_cg)
File(in_dir + "m.xml") >> m

# Newton solver params
prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-5
prm['newton_solver']['absolute_tolerance'] = 1e-5
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

pcs['k'] = 5e-3

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm
model_inputs['constants'] = pcs

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# Seconds per month
spm = pcs['spm']
# End time
T = 100.0 * spds
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
    model.write_xmls(['h', 'pfo', 'k', 'u_b'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
