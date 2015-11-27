
"""
Starts from steady state generated by sim_flat_high_melt_steady.py. Reduces melt,
conductivity, and sliding speed. Low melt variant. Conductivity lags one month
behind melt shutoff.
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
out_dir = "out_flat_bed_low_month_lag/"
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
File(in_dir + "h_low_steady.xml") >> h_init

# Sliding speed
u_b = Function(V_cg)
File(in_dir + "u_b.xml") >> u_b

# Seconds per year
spy = pcs['spy']
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

# Minimum conductivity
k_min = 5e-5
# Scaling parameter
a = 31220.5955507
# Parameter that controls lag of conductivity behind melt
b = 30.0 * spd
# Function that scales k proportionally to m
def k_scale(t):
  return a * m_scale(t - b)
  
# Create a function that scale u_b down to a maximum of 80 (m/a) in winter
c = (80.0 / spy) / u_b.vector().max()
def u_b_scale(t):
  return -m_scale(t) * (c - 1.0) + c

# Distributed melt
m = Function(V_cg)
File(in_dir + "m.xml") >> m
# Cut melt in half
m = project(0.5 * m, V_cg)

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
model_inputs['constants'] = pcs

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# Seconds per month
spm = pcs['spm']
# End time
T = 10.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Iteration count
i = 0

while model.t < T:  
  
  # Update the melt
  model.set_m(project(Constant(m_scale(model.t)) * m, V_cg))
  # Update the conductivity
  model.set_k(project(Constant(k_scale(model.t)) * m + Constant(k_min), V_cg))
  # Update the sliding speed
  model.set_u_b(project(Constant(u_b_scale(model.t)) * u_b, V_cg))  
  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 3 == 0:
    model.write_pvds(['h', 'u_b', 'm', 'pfo', 'k'])
    
  if i % 3 == 0:
    model.write_xmls(['h', 'pfo', 'k', 'u_b'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
