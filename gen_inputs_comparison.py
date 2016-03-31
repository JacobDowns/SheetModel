"""
This script generates all fields such as ice thickness, bed elevation, and melt
needed to run the model for an ice sheet on a flat bed. Compares directly to 
the first synthetic run in Mauro's paper.
"""

from dolfin import *
import numpy as np
from constants import *
from cr_tools import *

# Directory to write model inputs
out_dir = "inputs_channel/"
mesh = Mesh(out_dir + "mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)
V_cr = FunctionSpace(mesh, "CR", 1)

# Write inputs to a hdf5 file
f = HDF5File(mesh.mpi_comm(), out_dir + "inputs_channel_comparison.hdf5", 'w')
# Write the mesh to a file
f.write(mesh, "mesh")        
        
# Sliding speed
u_b = project(Constant(1e-6), V_cg)
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
    value[0] = sqrt((x[0] + 325.0) * h_max**2 / length)
    
# Surface
S = project(Surface(), V_cg)
# Bed elevation
B = project(Bed(), V_cg)
# Ice thickness
H = project(S - B, V_cg)


f.write(B, "B")
f.write(H, "H")

# Melt
day = 60.0 * 60.0 * 24.0        
class Melt(Expression):
  def eval(self,value,x):
    value[0] = max((0.14 - sqrt(x[0] * h_max**2 / length) * 1e-4) / day, 0.0)

m = project(Melt(), V_cg)
plot(m, interactive = True)
f.write(m, "m_0")

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

