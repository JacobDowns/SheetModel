"""
Starts from steady state generated by sim_trough_high_ref_steady.py and reduces
melt.
"""

import sys
from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '../')
from sheet_model import *
from constants import *

# Model input directory
in_dir = "inputs_high_melt/"
# Output directory
out_dir = "out_trough_high_ref_time/"
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
File(in_dir + "h_high_ref_steady.xml") >> h_init

# Melt input
m = Function(V_cg)
File(in_dir + "m.xml") >> m

# Seconds per month
spm = pcs['spm']
# Seconds per day
spd = pcs['spd']
# How long it takes for the melt to shut off completely
shutoff_length = 30.0 * spd
 
# Function that reduces melt to 0 over time
def m_scale(t):
  if t < 0.0:
    return 1.0
  if t <= shutoff_length:
    return cos((pi / (2.0 * shutoff_length)) * t)
  return 0.0

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

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# End time
T = 10.0 * spm
# Time step
dt = 60.0 * 60.0 * 12.0
# Iteration count
i = 0

while model.t < T:  
  
  # Update the melt
  model.set_m(project(Constant(m_scale(model.t)) * m, V_cg))
  
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

