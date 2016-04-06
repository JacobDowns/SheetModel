"""
This script generates all fields such as ice thickness, bed elevation, and melt
needed to run the model for an ice sheet overlying a trough. Low melt variant.
"""

from dolfin import *
import numpy as np
from constants import *
from cr_tools import *

# Directory to write model inputs
out_dir = "inputs_channel/inputs/"
mesh = Mesh("inputs_channel/mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)
V_cr = FunctionSpace(mesh, "CR", 1)

# Write inputs to a hdf5 file
f = HDF5File(mesh.mpi_comm(), out_dir + "inputs_trough_low.hdf5", 'w')
# Write the mesh to a file
f.write(mesh, "mesh")

# Melt
m = project(Expression("0.5 * (1.0 + (4.0 * (60000.0 - x[0]) / 60000.0)) / 31536000.0"), V_cg)
f.write(m, "m_0")

# Sliding speed
u_b = project(Expression("(50.0 + 250.0 * (60000.0 - x[0]) / 60000.0) / 31536000.0"), V_cg)
f.write(u_b, "u_b_0")


### Bed and surface functions

# Maximum ice thickness
h_max = 1500.
# Length of ice sheet 
length = 60e3
# Center of trough 
center = 10e3
# Maximum trough depth
depth = 500.0

def shape(x, s):
  return np.arctan(x / s) / (pi / 2.0)
  
def shape1(x, c, s):
  return (np.arctan((x - c) / s) / pi) + 0.5
      
class Bed(Expression):
  def eval(self,value,x):
    if x[0] < 5000.0:
      value[0] = 0.0
    else :
      X = x[0] - 5000.0
      Y = x[1]
      value[0] = -depth * e**(-((Y - center) / (2000.0 * shape(X, 500.0)) )**2.0) * shape(X, 7000.0)

class Surface(Expression):
  def eval(self,value,x):
    value[0] = sqrt((x[0] + 250.0) * h_max**2 / length)

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

f.write(boundaries, "boundaries")

# Initial sheet height
h = Function(V_cg)
h.interpolate(Constant(0.05))
f.write(h, "h_0")

# Initial channel height
S_0 = Function(V_cr)
f.write(S_0, "S_0")

# Compute the lengths of edges
edge_lens = Function(V_cr)
v_cr = TestFunction(V_cr)
L = assemble(v_cr('-') * dS + v_cr * ds).array()
edge_lens.vector()[:] = L
f.write(edge_lens, "edge_lens")

# Create a mask that's 0 on the boundary and 1 on the interior
M = assemble(v_cr * ds).array()
mask = Function(V_cr)
M[M != 0.0] = 2.0
M[M == 0.0] = 1.0
M[M == 2.0] = 0.0
mask.vector()[:] = M
#cr_tools.plot_cr(mask)
f.write(mask, "mask")

