from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world

# Model input directory
in_dir = "inputs_slope/"
# Output directory
out_dir = "out_comp_opt/"
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
#File("h_2299.xml") >> h_init

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Melt
u0 = Expression('sin(t)', t = 0)

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 0.96
prm['newton_solver']['relative_tolerance'] = 2e-3
prm['newton_solver']['absolute_tolerance'] = 1e-3
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 25

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['d_bcs'] = [bc]
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Time step
dt = 60.0 * 60.0

model.step_opt(dt)
model.write_pvds()

"""
model.phi_solver.solve_pde()
model.phi_solver.phi_apply_bounds()
model.update_phi()
model.write_pvds()"""
  

