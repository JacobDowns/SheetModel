"""
Model inputs for intercomparison test E1. 
"""

from dolfin import *
import numpy as np
from constants import *
from valley_mesh_helper import *

# Directory to write model inputs
in_dir = "../../inputs/"
out_dir = "../../inputs/E1/"
mesh = Mesh(in_dir + "mesh/mesh_valley.xml")


V_cg = FunctionSpace(mesh, "CG", 1)
V_tr = FunctionSpace(mesh, FiniteElement("Discontinuous Lagrange Trace", "triangle", 0))

# Write inputs to a hdf5 file
f = HDF5File(mesh.mpi_comm(), out_dir + "inputs_E1.hdf5", 'w')
# Write the mesh to a file
f.write(mesh, "mesh")

# Melt
m = interpolate(Constant(1.158e-6), V_cg)
f.write(m, "m_0")

# Sliding speed
u_b = interpolate(Constant(1e-6), V_cg)
f.write(u_b, "u_b_0")

k = interpolate(Constant(5e-3), V_cg)
f.write(k, "k_0")

### Bed and surface functions

      
class Bed(Expression):
  def eval(self,value,x):
    value[0] = bed(x[0], x[1], para_bench)

class Surface(Expression):
  def eval(self,value,x):
    value[0] = surface(x[0], x[1]) + 1.0

# Surface
S = project(Surface(), V_cg)
# Bed elevation
B = project(Bed(), V_cg)
# Ice thickness
H = project(S - B, V_cg)

f.write(B, "B")
f.write(H, "H")

# Initial pressure
# Potential at 0 pressure
phi_m = project(pcs['rho_w'] * pcs['g'] * B, V_cg)
# Ice overburden pressure
p_i = project(pcs['rho_i'] * pcs['g'] * H, V_cg)
# Potential at overburden pressure
phi_0 = project(phi_m + p_i, V_cg)
f.write(phi_0, "phi_0")

### Create a facet function with marked boundaries

# Margin
class MarginSub(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], 0.0)
    
ms = MarginSub()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
ms.mark(boundaries, 1)

plot(boundaries, interactive = True)

f.write(boundaries, "boundaries")

# Initial sheet height
h = Function(V_cg)
h.interpolate(Constant(0.06))
f.write(h, "h_0")

# Initial channel height
S_0 = Function(V_tr)
f.write(S_0, "S_0")

