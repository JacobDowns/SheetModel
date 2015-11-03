from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world
from numpy import *

""" Generates a steady state for diurnal simulations."""

# Model input directory
in_dir = "inputs_sliding_law/"
# Output directory
out_dir = "out_diurnals_steady/"
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
File(in_dir + "h_steady.xml") >> h_init

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Driving stress
tau_b = Function(V_cg)
File(in_dir + "tau_b.xml") >> tau_b

# Load point source melt function 
m = Function(V_cg) 
File(in_dir + "m_point2.xml") >> m

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

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
model_inputs['m'] = m
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
T = 200.0 * spd
# Time step
dt = 60.0 * 60.0 * 4.0
# Iteration count
i = 0
# Seconds per year
spy = pcs['spy']
# Somewhat arbitrary constant for sliding law
C = Constant(0.1e-8 / spy)
u_b_reg = Constant(100.0)

while model.t < T:
  if MPI_rank == 0: 
    current_time = model.t / spd
    print "Current Time: " + str(current_time)
  
  model.step(dt)
  
  if i % 3 == 0:
    model.write_pvds(['h', 'pfo', 'u_b'])
    
  if i % 6 == 0:
    model.write_xmls(['h', 'u_b', 'phi'])
  
  if MPI_rank == 0: 
    print
    
  # Update the sliding velocity
  model.update_u_b(project(C * (tau_b**3 / (model.N + u_b_reg)), V_cg))
  
  i += 1