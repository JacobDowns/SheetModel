"""
Model inputs for intercomparison test E1. 
"""

from dolfin import *
import numpy as np
from constants import *
from valley import *

# Bedrock parameters
params = [0.05, 0.0, -0.1, -0.5, -0.7]

# Load the valley outline mesh
mesh = Mesh("../../inputs/mesh/mesh_valley.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Melt
m = interpolate(Constant(1.158e-6), V_cg)
# Sliding speed
u_b = interpolate(Constant(1e-6), V_cg)
# Conductivity
k = interpolate(Constant(5e-3), V_cg)

h = interpolate(Constant(0.01), V_cg)

# Margin boundary 

class MarginSub(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], 0.0)
    
ms = MarginSub()
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
ms.mark(boundaries, 1)

## Write input files
for n in range(5):
  print n
  
  out_file = "../../inputs/E/input_E" + str(n + 1) + ".hdf5"
  
  bed_param = params[n]  
  f = HDF5File(mesh.mpi_comm(), out_file, 'w')
  f.write(mesh, "mesh")
  f.write(m, "m_0")
  f.write(u_b, "u_b_0")
  f.write(k, "k_0")  
  f.write(boundaries, "boundaries")
  f.write(h, "h_0")
  
  ## Bed and surface functions
  class Bed(Expression):
    def eval(self,value,x):
      bed, thickness = valley(np.array([x[0]]), np.array([x[1]]), bed_param)
      value[0] = bed[0]
  
  class Thickness(Expression):
    def eval(self,value,x):
      bed, thickness = valley(np.array([x[0]]), np.array([x[1]]), bed_param)
      value[0] = thickness[0]
  
  # Surface
  H = project(Thickness(degree = 1), V_cg)  
  f.write(H, "H")
  # Bed elevation
  B = project(Bed(degree = 1), V_cg)
  f.write(B, "B")
  
  File('B' + str(n) + '.pvd') << B
  File('H' + str(n) + '.pvd') << H
  
  # Potential at 0 pressure
  phi_m = project(pcs['rho_w'] * pcs['g'] * B, V_cg)
  # Ice overburden pressure
  p_i = project(pcs['rho_i'] * pcs['g'] * H, V_cg)
  # Potential at overburden pressure
  phi_0 = project(phi_m + p_i, V_cg)
  f.write(phi_0, "phi_0")
