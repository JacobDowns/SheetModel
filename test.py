from dolfin import *
from constants import *
from dolfin import MPI, mpi_comm_world
from numpy import *

# Model input directory
in_dir = "inputs_is/"
# Output directory
out_dir = "out_is/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())

"""
# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

run_name = "test"
f = HDF5File(mesh.mpi_comm(), run_name + ".hdf5", 'w')
f.write(mesh, "mesh")

h = Function(V_cg)
out = File("test/h.xdmf")

h.interpolate(Constant(1.0))
f.write(h, "h")
#out << h


h.interpolate(Constant(2.0))
f.write(h, "h", 2.0)
out << h

h.interpolate(Constant(3.0))
f.write(h, "h", 2.0)
out << h
#f.write(array([1.0]), "ts")
#f.write(array([2.0]), "ts")"""




f = HDF5File(mpi_comm_world(), "test.hdf5", 'r')

# Load mesh
mesh = Mesh()
f.read(mesh, "mesh", False)


#plot(mesh, interactive = True)

count = f.attributes("h")
print count 

quit()
V_cg = FunctionSpace(mesh, "CG", 1)
h = Function(V_cg)
f.read(h, "h/vector_1")

#plot(h, interactive = True) 

attr = f.attributes("h/vector_" + str(count - 1))


