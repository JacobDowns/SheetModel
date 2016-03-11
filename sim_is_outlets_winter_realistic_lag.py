"""
Simulates winter pressure on IS. Conductivity and sliding vary over time.
Outlets bcs. Lag time is 1 day.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_is_outlets_winter_realistic_lag'
# Input file
input_file = 'inputs/steady_realistic/outlets_steady_realistic.hdf5'
# Load the input file
inputs = HDF5File(mpi_comm_world(), input_file, 'r')


### Create outlet boundary conditions
mesh = Mesh()
inputs.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1)

# Load potential at 0 pressure
B = Function(V_cg)
inputs.read(B, "B")

def outlet_boundary(x, on_boundary):
  # These two outlet points are based on Google Earth -- places where it looks
  # like water is flowing out
  out1_x = -491254.66
  out1_y = -2461159.0
  out2_x = -491839.3
  out2_y = -2472998.0
  
  cond1 = (abs(x[0] - out1_x) < 150.0) and (abs(x[1] - out1_y) < 150.0)
  cond2 = (abs(x[0] - out2_x) < 150.0) and (abs(x[1] - out2_y) < 150.0)
  return cond1 or cond2

bc = DirichletBC(V_cg, pcs['rho_w'] * pcs['g'] * B, outlet_boundary, "pointwise")


### Initialize model

# Use a smaller conductivity
pcs['k'] = 7e-3
model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['d_bcs'] = [bc]
model_inputs['constants'] = pcs
model_inputs['opt_params'] = {'tol' : 1e-2, 'scale' : 30}
model = SheetModel(model_inputs)

m = Function(model.V_cg)
m.assign(model.m)
u_b = Function(model.V_cg)
u_b.assign(model.u_b)


### Set up scaling functions for time dependent m, k, and u_b

# Seconds per day
spd = pcs['spd']
# Seconds per year
spy = pcs['spy']
# How long it takes for the melt to shut off completely
shutoff_length = 30.0 * spd
 
# Function that reduces melt to 0 over time
def m_scale(t):
  if t < 0.0:
    return 1.0
  if t <= shutoff_length:
    return cos((pi / (2.0 * shutoff_length)) * t)
  return 0.0

# Maximum conductivity
k_max = model.pcs['k']
# Minimum conductivity
k_min = 7e-5
# Scaling parameter that sets the maximum possible conductivity
a = (k_max - k_min) / m.vector().max()
# Lag time
b = 7.0*spd
# Function that scales k proportionally to m
def k_scale(t):
  return a * m_scale(t-b)
  
# Create a function that scales u_b down to a maximum of 80 (m/a) in winter
c = (100.0 / pcs['spy']) / u_b.vector().max()
def u_b_scale(t):
  return -m_scale(t) * (c - 1.0) + c


### Run the simulation

# Seconds per day
spd = pcs['spd']
# Seconds per month
spm = pcs['spm']
# End time
T = 7.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Iteration count
i = 0

while model.t < T:  
  
  # Update the melt
  model.set_m(project(Constant(m_scale(model.t)) * m, model.V_cg))
  # Update the conductivity
  model.set_k(project(Constant(k_scale(model.t)) * m + Constant(k_min), model.V_cg))
  # Update the sliding speed
  model.set_u_b(project(Constant(u_b_scale(model.t)) * u_b, model.V_cg))  
  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 3 == 0:
    model.write_pvds(['h', 'pfo'])
    
  if i % 3 == 0:
    model.checkpoint(['pfo', 'h', 'k', 'm', 'u_b'])
  
  if MPI_rank == 0: 
    print
    
  i += 1
