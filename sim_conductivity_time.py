from dolfin import *
from dolfin_adjoint import *
from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world
from pylab import *

""" 
This simulation starts from the steady state generated by sim_sliding_law_time.py.
During the winter hydraulic condictivity is lowered from it's default value."""

# Model input directory
in_dir = "inputs_sliding_law/"
# Output directory
out_dir = "out_conductivity_time/"
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
File(in_dir + "h_sliding_steady.xml") >> h_init

# Load the boundary facet function
boundaries = FacetFunction('size_t', mesh)
File(in_dir + "boundaries.xml") >> boundaries

# Load potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m

# Driving stress
tau_b = Function(V_cg)
File(in_dir + "tau_b.xml") >> tau_b

# Function that modifies the melt function over time
# Seconds per year
spy = pcs['spy']
# Parameters for logistic functions
t0 = 0.825 * spy
t1 = 0.175 * spy
# Value to scale the logistic curve
M = 1.0 / (1.0 + e**(-7e-7 * t1))

def m_time(t):
  melt1 = 1.0 / (1.0 + e**(-7e-7 * (t - t0)))
  melt2 = 1.0 / (1.0 + e**(-7e-7 * (-t + t1)))  
  return max(melt1, melt2) * (1.0 / M)

# Define a time varying hydraulic conductivity
# Maximum and minimum conductivities
min_k = 1e-3
max_k = 5e-3

def k_time(t):
  k1 = 1.0 / (1.0 + e**(-7e-7 * (t - t0)))
  k2 = 1.0 / (1.0 + e**(-7e-7 * (-t + t1))) 
  kval = max(k1, k2) * ((max_k - min_k) / M) + min_k
  return kval   

  
# Load initial melt function
m = Function(V_cg) 
File(in_dir + "m.xml") >> m
# Copy melt function
m1 = Function(V_cg)
m1.assign(m)

# Initial sliding velocity
u_b = Function(V_cg)
File(in_dir + "u_b_steady.xml") >> u_b

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
model_inputs['u_b'] = u_b
model_inputs['k'] = k
model_inputs['d_bcs'] = [bc]
model_inputs['out_dir'] = out_dir
model_inputs['newton_params'] = prm
model_inputs['constants'] = pcs
model_inputs['opt_params'] = {'tol' : 0.5e-7, 'scale' : 15.0}

# Create the sheet model
model = SheetModel(model_inputs, in_dir)


### Run the simulation

# Seconds per year
spy = pcs['spy']
# Seconds per day
spd = pcs['spd']
# End time
T = 1.5 * spy
# Time step
dt = 60.0 * 60.0 * 4.0
# Iteration count
i = 0
# Somewhat arbitrary constant for sliding law
C = Constant(0.33e-8 / spy)

while model.t < T:
  if MPI_rank == 0: 
    current_time = model.t / spd
    print "Current Time: " + str(current_time)
  
  model.step(dt)
  
  if i % 3 == 0:
    model.write_pvds(['h', 'pfo', 'u_b', 'm', 'k'])
    
  if i % 6 == 0:
    model.write_xmls(['pfo', 'h', 'u_b'])
  
  if MPI_rank == 0: 
    print
    
  # Update the sliding velocity
  model.update_u_b(project(Constant(C) * (tau_b**3 / model.N), V_cg))
  # Update the melt
  model.update_m(project(Constant(m_time(model.t)) * m1, V_cg))
  # Update conductivity
  model.update_k(project(Constant(k_time(model.t)), V_cg))
  
  i += 1