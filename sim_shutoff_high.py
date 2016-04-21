"""
Gammut of high melt shutoff time tests.
"""

from dolfin import *
from constants import *
from sheet_model import *
from dolfin import MPI, mpi_comm_world
from scale_functions import *

MPI_rank = MPI.rank(mpi_comm_world())


### Load the input file 

input_file = 'inputs_sheet/steady/real_high_steady.hdf5'
k_min = 5e-5
k_max = 5e-3
scale_functions = ScaleFunctions(input_file, k_min, k_max, u_b_max = 80.0)


### Setup the model

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 30

# Seconds per day
spd = pcs['spd']
model_inputs = {}
model_inputs['constants'] = pcs
model_inputs['input_file'] = input_file
model_inputs['newton_params'] = prm

shutoff_times = [spd, 7.0 * spd, 30.0 * spd]
names = ['day', 'week', 'month']


### Run the simulation

def run_sim(name):

  # Seconds per month
  spm = pcs['spm']
  # End time
  T = 8.0 * spm
  # Time step
  dt = 60.0 * 60.0 * 8.0
  # Iteration count
  i = 0
  
  while model.t < T:  
    if name == 'day' and model.t <= spd:
      dt = 60.0 * 60.0 * 1.0
    elif name == 'week' and model.t <= 7.0 * spd:
      dt = 60.0 * 60.0 * 1.0
    else :
      dt = 60.0 * 60.0 * 8.0
      
    # Update the melt
    model.set_m(scale_functions.get_m(model.t))
    # Update the conductivity
    model.set_k(scale_functions.get_k(model.t))
    # Update the sliding speed
    model.set_u_b(scale_functions.get_u_b(model.t))  
    
    if MPI_rank == 0: 
      current_time = model.t / spd
      print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
    
    model.step(dt)
    
    if i % 1 == 0:
      model.write_pvds(['pfo', 'h'])
      
    if i % 1 == 0:
      model.checkpoint(['h', 'pfo', 'k', 'm', 'u_b'])
    
    if MPI_rank == 0: 
      print
      
    i += 1
    
    
### Run simulations with different lag times 

i = 0
for to in shutoff_times:
  scale_functions.shutoff_length = to
  name = names[i]
  model_inputs['out_dir'] = 'paper_results/real_shutoff_high/' +  name
  model_inputs['checkpoint_file'] = 'shut_high_' + name
  model = SheetModel(model_inputs)
  run_sim(name)
  i += 1
