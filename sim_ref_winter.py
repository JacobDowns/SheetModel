"""
Reference simulation over winter.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 30

model_inputs = {}
pcs['k'] = 2e-2
pcs['alpha'] = 2
model_inputs['input_file'] = 'inputs/steady_ref/steady_ref.hdf5'
model_inputs['out_dir'] = 'out_ref_winter/'
model_inputs['constants'] = pcs
model_inputs['newton_params'] = prm

# Create the sheet model
model = SheetModel(model_inputs)


### Melt scale function 

# Copy the summer high melt function
m = Function(model.V_cg)
m.assign(model.m)

# Seconds per day
spd = pcs['spd']
# How long it takes for the melt to shut off completely
shutoff_length = 30.0 * spd
 
# Function that reduces melt to 0 over time
def m_scale(t):
  if t < 0.0:
    return 1.0
  if t <= shutoff_length:
    return cos((pi / (2.0 * shutoff_length)) * t)
  return 0.0


### Run the simulation

# Seconds per month
spm = pcs['spm']
# Seconds per day
spd = pcs['spd']
# End time
T = 4.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Iteration count
i = 0

while model.t < T:  
  # Update the melt
  model.set_m(project(Constant(m_scale(model.t)) * m, model.V_cg))
  
  if MPI_rank == 0: 
    current_time = model.t / spd
    #print "Current Time: " + str(current_time)
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  
  if i % 1 == 0:
    model.write_pvds(['pfo', 'h'])
    
  if i % 1 == 0:
    model.checkpoint(['m','pfo'])
  
  if MPI_rank == 0: 
    print
    
  i += 1