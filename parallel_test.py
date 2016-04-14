# -*- coding: utf-8 -*-
"""
A test of the PDE solve in parallel.
"""

from dolfin import *
from constants import *
from channel_model import *
from dolfin import MPI, mpi_comm_world
from scale_functions import *

MPI_rank = MPI.rank(mpi_comm_world())


### Load the input file 

input_file = 'inputs_channel/steady/real_trough_steady_high.hdf5'
k_min = 5e-5
k_max = 5e-3
scale_functions = ScaleFunctions(input_file, k_min, k_max, u_b_max = 80.0)


### Setup the model

prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-12
prm['newton_solver']['absolute_tolerance'] = 1e-12
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 50


# Seconds per day
spd = pcs['spd']
model_inputs = {}
pcs['k_c'] = 1e-1
model_inputs['constants'] = pcs
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = 'parallel_test'
model_inputs['newton_params'] = prm
model = ChannelModel(model_inputs)

# Update the melt
model.set_m(scale_functions.get_m(model.t))
# Update the conductivity
model.set_k(scale_functions.get_k(model.t))
# Update the sliding speed
model.set_u_b(scale_functions.get_u_b(model.t))  


### Do one pressure solve

# Time step
dt = 60.0 * 60.0 * 8.0
model.step(dt)

File("parallel_test/phie.pvd") << model.phi
