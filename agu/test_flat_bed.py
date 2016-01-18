import sys
from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '../')
from sheet_model import *
from constants import *
from pylab import *

# Model input directory
in_dir = "inputs_flat_bed/"
# Output directory
out_dir = "out_flat_bed/"
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
File(in_dir + ")

# Seconds per day
spd = pcs['spd']
# How long it takes for the melt to shut off completely
shutoff_length = 30.0 * spd
 
# Function that reduces melt to 0 over time
def m_scale(t):
  if t <= shutoff_length:
    return cos((pi / (2.0 * shutoff_length)) * t)
  return 0.0

# Minimum conductivity
k_min = 5e-5
# Maximum conductivity
k_max = 5e-3
# Scaling parameter
a = (k_max - k_min) / m.vector().max()
# Parameter that controls lag of conductivity behind melt
b = 0.0

def k_scale(t):
  return a * m_scale(t - b)
  
# Create a function that scale u_b down to a maximum of 80 (m/a) in winter
c = (80.0 / spy) / u_b.vector().max()
def u_b_scale(t):
  return -m_scale(t) * (c - 1.0) + c

t0 = 2.0 * spm

# Distributed melt function (low melt)
m = Function(V_cg)
File(in_dir + "m.xml") >> m

# Moulin melt function
m_moulin = Function(V_cg)
File(in_dir + "m_moulin.xml") >> m_moulin

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
# Time step
dt = 60.0 * 60.0 * 8.0

model.step(dt)
plot(model.pfo, interactive = True)



