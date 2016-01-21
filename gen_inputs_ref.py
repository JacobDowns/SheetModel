"""
This script generates all fields such as ice thickness, bed elevation, and melt
needed to run the model for an ice sheet on a flat bed. High melt variant.
"""

from dolfin import *
import numpy as np
from constants import *

# Directory to write model inputs
out_dir = "inputs/"
mesh = Mesh(out_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Write inputs to a hdf5 file
f = HDF5File(mesh.mpi_comm(), out_dir + "inputs_ref.hdf5", 'w')
# Write the mesh to a file
f.write(mesh, "mesh")

# Melt
m = project(Expression("(1.0 + (4.0 * (60000.0 - x[0]) / 60000.0)) / 31536000.0"), V_cg)
f.write(m, "m_0")

# Sliding speed
u_b = project(Expression("(50.0 + 250.0 * (60000.0 - x[0]) / 60000.0) / 31536000.0"), V_cg)
f.write(u_b, "u_b_0")


### Bed and surface functions

# Maximum ice thickness
h_max = 1500.
# Length of ice sheet 
length = 60e3
      
class Bed(Expression):
  def eval(self,value,x):
    value[0] = 0.0

class Surface(Expression):
  def eval(self,value,x):
    value[0] = sqrt((x[0] + 50.0) * h_max**2 / length)

# Surface
S = project(Surface(), V_cg)
# Bed elevation
B = project(Bed(), V_cg)
# Ice thickness
H = project(S - B, V_cg)

f.write(B, "B")
f.write(H, "H")


### Create a facet function with marked boundaries

# Margin
class MarginSub(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], 0.0)
    
ms = MarginSub()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
ms.mark(boundaries, 1)

f.write(boundaries, "boundaries")

# Initial sheet height
h = Function(V_cg)
h.interpolate(Constant(0.05))
f.write(h, "h_0")
