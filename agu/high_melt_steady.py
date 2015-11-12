"""
Runs the trough to a steady state (runs with constant melt input for 3 months),
then reduces melt input. Conductivity is scaled in the same way as melt.
"""

import sys
from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '/home/jake/schoof/')
from sheet_model import *
from constants import *
from pylab import *

# Model input directory
in_dir = "inputs_high_melt/"
# Output directory
out_dir = "out_high_melt_steady/"
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

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Melt function
m = Function(V_cg)
File(in_dir + "m.xml") >> m

spy = pcs['spy']
spm = spy / 12.0

# Function that reduces melt to 0 over time
def m_scale(t):
  # Keep the melt steady for 3 months
  if t <= 3.0 * spm:
    return 1.0
    
  # After 3 months decrease melt
  t = t - 3.0 * spm  
  
  return 1.0 / (1.0 + e**(1e-6 * (t - 3.0 * spm)))
  
# Function that reduces k proportionally to m
def k_scale(t):
  # Maximum k value
  k_max = 5e-3
  # Minimum k value
  k_min = 5e-4
  # Lag of conductivity behind melt
  delay = 0.0
  
  return (k_max - k_min) * m_scale(t - delay) + k_min

# Function that scale the sliding speed
u_b = Function(V_cg)
File(in_dir + "u_b.xml") >> u_b
dolfin.plot(u_b * spy, interactive = True)
quit()

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

# Use a slightly lower conductivity than the default
pcs['k'] = 1e-3

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-3
prm['newton_solver']['absolute_tolerance'] = 1e-3
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['d_bcs'] = [bc]
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm
model_inputs['m'] = m
model_inputs['constants'] = pcs

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
    print "Current Time: " + str(current_time)
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds()
    
  if i % 1 == 0:
    model.write_xmls()
  
  if MPI_rank == 0: 
    print
    
  i += 1

