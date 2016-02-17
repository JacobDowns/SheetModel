"""
Solves for pressure with the PDE then the optimization procedure for comparison.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Output directory
out_dir = 'out_is_outlets_winter'
# Input file
input_file = 'inputs/inputs_is_outlets_steady.hdf5'
# Load the input file
inputs = HDF5File(mpi_comm_world(), input_file, 'r')

### Create outlet boundary conditions
mesh = Mesh()
inputs.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1)

phi_min = project(Constant(-100000000.0), V_cg)
phi_max = project(Constant(100000000.0), V_cg)

### Initialize model

# Use a smaller conductivity
pcs['k'] = 5e-3

model_inputs = {}
model_inputs['input_file'] = input_file
model_inputs['out_dir'] = out_dir
model_inputs['checkpoint_file'] = 'outlets_is_winter'
model_inputs['constants'] = pcs
model_inputs['phi_min'] = phi_min
model_inputs['phi_max'] = phi_max


# Create the sheet model
model = SheetModel(model_inputs)

# Solve for phi using the pde
#model.phi_solver.solve_pde()

# Solve with optimization procedure
model.phi.assign(model.phi_m)
phi1 = Function(V_cg)
model.phi_solver.solve_opt()
phi1.assign(model.phi)

# Solve with pde
model.phi.assign(model.phi_m)
phi2 = Function(V_cg)
model.phi_solver.solve_pde()
phi2.assign(model.phi)

File("opt_test.pvd") << project(abs(phi1 - phi2) / abs(phi1), V_cg)