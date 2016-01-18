from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world
import numpy as np

"""
Starts from the steady state from sim_sliding_steady.py. This is a trough
simulation where everything is fixed except the bump height which is varied.
"""

# Model input directory
in_dir = "inputs_sliding_law/"
# Output directory
out_dir = "out_sliding_sensitivity1/"
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
h = Function(V_cg)
h.assign(h_init)

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Enforce 0 pressure bc at margin
bc = DirichletBC(V_cg, phi_m, boundaries, 1)

# Sliding speed
u_b = Function(V_cg)
spy = pcs['spy']
u_b.interpolate(Constant(100.0 / spy))


pcs['h_r'] = 0.075
pcs['k'] = 5e-3

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

# Time step
dt = 60.0 * 60.0 * 4.0

hrs = np.linspace(0.05, 0.5, 40)

t_10 = []
t_20 = []
t_50 = []
o_10 = []
o_20 = []
o_50 = []



model.step(dt)
plot(model.pfo, interactive = True)
"""

for hr in hrs:
  model.pcs['h_r']
  model.step(dt)
  
  #plot(model.pfo, interactive = True)
  
  # Reset the solution
  model.h.assign(h)
  
  model.write_xmls(['pfo', 'phi'])"""
  
  
