from dolfin import *
from dolfin_adjoint import *
#from sheet_model import *
from constants import *
from dolfin import MPI, mpi_comm_world


# Model input directory
in_dir = "inputs_is/"
# Output directory
out_dir = "out_is/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())


# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

run_name = "test"
f = HDF5File(mesh.mpi_comm(), run_name + ".hdf5", 'w')
#f.write(mesh, "mesh")


h = Function(V_cg)

h.interpolate(Constant(1.0))
f.write(h, "h", 1.0)

h.interpolate(Constant(2.0))
f.write(h, "h", 2.0)

a = "a" + 1

h.interpolate(Constant(3.0))
f.write(h, "h", 3.0)

h.interpolate(Constant(4.0))

attrs = f.attributes("h")



"""
f = HDF5File(mpi_comm_world(), "test.hdf5", 'r')

# Load mesh
mesh = Mesh()
f.read(mesh, "mesh", False)


# Load hs
try :
  print f.attributes("h")
except:
  print "sfas"

print dir(f)

quit()



#plot(mesh, interactive = True)

count = f.attributes("h")['count']

V_cg = FunctionSpace(mesh, "CG", 1)
h = Function(V_cg)
f.read(h, "h/vector_" + str(count - 1))

plot(h, interactive = True) 

attr = f.attributes("h/vector_" + str(count - 1))"""


