"""
IS steady state simulation.
"""

import sys
from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
sys.path.insert(0, '../')
from sheet_model import *
from constants import *


# Model input directory
in_dir = "inputs_is/"
# Output directory
out_dir = "out_test/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())
# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

g = Function(V_cg, name = "g")

mesh_name = "mesh"
f = HDF5File(mesh.mpi_comm(), out_dir + mesh_name + ".hdf5", 'w')

f.write(mesh, mesh_name)
f.write(g, "waffles/g")

g_out = File(out_dir + "g.xdmf")
g_out << (g, 1.0)
g_out << (g, 2.0)
g_out << (g, 3.0)

"""
mesh = Mesh()
f = HDF5File(mpi_comm_world(), out_dir + "mesh.hdf5", 'r')
f.read(mesh, "mesh", False)

V_cg = FunctionSpace(mesh, "CG", 1)

plot(mesh, interactive = True)

g = Function(V_cg)
f.read(g, "mesh/g")

plot(g, interactive = True)"""