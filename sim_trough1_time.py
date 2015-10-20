from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world
from numpy import *

# Model input directory
in_dir = "inputs_trough1/"
# Output directory
out_dir = "out_trough1_time/"
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

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

# Use a slightly lower conductivity than the default
pcs['k'] = 5e-3

# Time dependent melt and sliding velocity
spy = pcs['spy']
# 5 - 1 meters summer to 0 meters winter
m = Expression('(1.0 + (4.0 * (60000.0 - x[0]) / 60000.0)) * 0.5 * (cos(((2.0 * pi) / spy) * t) + 1.0) / spy', pi = pi, spy = spy, t = 0.0)
# 200 - 50 meters summer to 150 - 0 meters winter
u_b = Expression("(25.0 * (cos(((2.0 * pi) / spy) * t) + 1.0) + 150.0 * (60000.0 - x[0]) / 60000.0) / spy", pi = pi, spy = spy, t = 0.0)
  
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
model_inputs['u_b'] = u_b
model_inputs['d_bcs'] = [bc]
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm
model_inputs['constants'] = pcs
model_inputs['opt_params'] = {'tol' : 1e-7, 'scale' : 10.0}

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 1.5 * spy
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
    model.write_pvds(['u_b', 'm', 'h', 'pfo'])
    
  if i % 6 == 0:
    model.write_xmls()
  
  if MPI_rank == 0: 
    print
    
  i += 1
