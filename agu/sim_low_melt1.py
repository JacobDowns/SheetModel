"""
Runs the trough to a steady state (runs with constant melt input for 3 months),
then reduces melt input. Conductivity is scaled in the same way as melt.
"""

import sys
from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '../')
from sheet_model import *
from constants import *
from pylab import *

# Model input directory
in_dir = "inputs_low_melt/"
# Output directory
out_dir = "out_low_melt_steady/"
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

# Distributed melt function (low melt)
m = Function(V_cg)
File(in_dir + "m.xml") >> m

# Moulin melt function
m_moulin = Function(V_cg)
File(in_dir + "m_moulin.xml") >> m_moulin
# Cut the moulin melt rate in half
m_moulin = project(m_moulin * 0.5, V_cg)

# Sliding speed
u_b = Function(V_cg)
File(in_dir + "u_b.xml") >> u_b 

spy = pcs['spy']
spm = pcs['spm']

# Function that reduces melt to 0 over time
def m_scale(t):
  # Keep the melt steady for 2 months
  if t <= 3.0 * spm:
    return 1.0
  # After 3 months decrease melt
  t = t - 3.0 * spm  
  
  return 1.0 / (1.0 + e**(2e-6 * (t - 2.0 * spm)))
  
# Create a function that scales k proportionally to m
# Minimum conductivity
k_min = 5e-5
# Scaling parameter
a = 31220.5955507
# Parameter that controls lag of conductivity behind melt
b = 0.0

def k_scale(t):
  return a * m_scale(t - b)
  
# Create a function that scale u_b down to a maximum of 80 (m/a) in winter
c = (80.0 / spy) / u_b.vector().max()
def u_b_scale(t):
  return -m_scale(t) * (c - 1.0) + c


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
model_inputs['m'] = m_moulin
model_inputs['constants'] = pcs

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 3.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Irataion count
i = 0

while model.t < T:
  # Update the melt
  model.set_m(project(Constant(m_scale(model.t)) * m_moulin, V_cg))
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

