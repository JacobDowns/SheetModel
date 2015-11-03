from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world

"""
Starts from the steady state from sim_sliding_steady.py. This is a trough
simulation where melt is fixed and sliding is varied with time but is uniform
in space.
"""

# Model input directory
in_dir = "inputs_sliding_law/"
# Output directory
out_dir = "out_sliding_test_time/"
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
File(in_dir + "h_sliding_test_steady.xml") >> h_init

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

# Seconds per month
spm = pcs['spy'] / 12.0
# Length of each transition
L = spm * 2.0

# Amount to scale the sliding speed
def u_b_scale(t):
  # Cycle e
  if t <= L:
    return 1.0 - (t / (L / 0.9))
  elif t <= 2.0 * L:
    return 0.1
  elif t <= 3.0 * L:
    return 0.1 + (t - 2.0 * L) / (L / 0.9)
  else:
    return 1.0

# Sliding speed
u_b = Function(V_cg)
spy = pcs['spy']
u_b.interpolate(Constant(100.0 / spy))
# Copy of sliding speed
u_b1 = Function(V_cg)

# Use a slightly lower conductivity than the default
pcs['k'] = 5e-3
# Bump height
pcs['h_r'] = 0.1

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-5
prm['newton_solver']['absolute_tolerance'] = 1e-5
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['u_b'] = u_b
model_inputs['d_bcs'] = [bc]
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm
model_inputs['constants'] = pcs
model_inputs['opt_params'] = {'tol' : 0.5e-7, 'scale' : 15.0}

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 8.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Iteration count
i = 0

while model.t < T:
  if MPI_rank == 0: 
    current_time = model.t / spd
    print "Current Time: " + str(current_time)
  
  model.step(dt)
  
  if i % 3 == 0:
    model.write_pvds(['h', 'pfo', 'u_b'])
    
  if i % 1 == 0:
    model.write_xmls()
  
  if MPI_rank == 0: 
    print
  
  model.update_u_b(project(Constant(100.0 * u_b_scale(model.t) / spy), V_cg))
  i += 1