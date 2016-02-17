from dolfin import *
from dolfin import MPI, mpi_comm_world
import numpy as np

MPI_rank = MPI.rank(mpi_comm_world())

mesh = Mesh('meshes/mesh_60_20.xml')

V_cg = FunctionSpace(mesh, 'CG', 1)

indexes = np.array(range(V_cg.dim()), dtype = np.intc)


f = Function(V_cg)
f.vector().set_local(np.array(range(len(f.vector().array())), dtype = np.float_))


v = TestFunction(V_cg)
F = f * v * dx

x = Vector()
assemble(F).gather(x, indexes)

print x.array()

"""

if MPI_rank == 0:
  #print MPI_rank
  f.vector().gather(x, indexes)
  print MPI_rank
  
print MPI_rank

quit()

global_indexes = V_cg.dofmap().tabulate_local_to_global_dofs()
print(len(global_indexes), len(f.vector().array()))


global_dofs = V_cg.dofmap().dofs()
f.vector().set_local(x[global_dofs])
f.vector().apply("insert")

print f.vector().array()"""


