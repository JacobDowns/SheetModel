from dolfin import *
from dolfin import MPI, mpi_comm_world
import numpy as np

MPI_rank = MPI.rank(mpi_comm_world())

mesh = Mesh('meshes/mesh_60_20.xml')

V_cg = FunctionSpace(mesh, 'CG', 1)

f = Function(V_cg)
v = TestFunction(V_cg)


F = f * v * dx



indexes = np.array(range(V_cg.dim()), dtype = np.intc)
x = Vector()

f.vector().gather(x, indexes)

print len(assemble(F).array())


print x.array()
print len(x)