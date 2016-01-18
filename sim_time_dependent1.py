from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world

# Model input directory
in_dir = "inputs_trough/"
# Output directory
out_dir = "out_time_dependent/"
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
# Load the steady state trough sheet height
File("out_trough_steady/checkpoint/h_224.xml") >> h_init

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Time dependent melt rate
# Seconds per year
spy = pcs['spy']
m = Expression('(2.0 + cos( ((2.0 * pi) / spy) * t)) / spy', pi = pi, spy = spy, t = 0.0)

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

# Use a slightly lower conductivity than the default
pcs['k'] = 5e-3

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-3
prm['newton_solver']['absolute_tolerance'] = 1e-3
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 30

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
T = 1.5 * spy
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
    model.write_pvds(set(['phi', 'pfo', 'h', 'm']))
    
  if i % 3 == 0:
    model.write_xmls(set(['phi', 'h']))
  
  if MPI_rank == 0: 
    print
    
  i += 1
