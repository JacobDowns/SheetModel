from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world

# Model input directory
in_dir = "inputs_slope/"
# Checkpoint directory
check_dir = "out/checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())

# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)


### Model inputs

# Initial sheet height
h_init = Function(V_cg)
File(check_dir + "h_2380.xml") >> h_init

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['d_bcs'] = [bc]
model_inputs['out_dir'] = "stuff"

# Create the sheet model
model = SheetModel(model_inputs, in_dir)

model.phi_solver.solve_opt(2e-7)
model.update_phi()

plot(model.phi, interactive = True)
File("comparison.pvd") << model.pfo
