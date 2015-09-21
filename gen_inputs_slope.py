"""
This script generates all fields such as ice thickness, bed elevation, and melt
needed to run the model for an ice sheet with a sloped bed.
"""

from dolfin import *
from parallel_map_gen import *
from constants import *

# Directory to write model inputs
out_dir = "inputs_slope/"
# Directory to write parallel maps
maps_dir = out_dir + "maps"

mesh = Mesh("inputs_slope/mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)
V_cr = FunctionSpace(mesh, "CR", 1)

### Create the parallel maps
# Create the parallel map generator
map_gen = ParallelMapGen(mesh)
# Write maps to a folder
map_gen.write_maps(maps_dir)


### Create a melt function (m/s)
m = project(Expression("(1.0 + 1.5 * (50000.0 - x[0]) / 50000.0) / 31536000.0"), V_cg)
File(out_dir + "m.xml") << m
File(out_dir + "m.pvd") << m


#### Create a sliding sliding speed (m/s)
u_b = project(Expression("(10.0 + 240.0 * (50000.0 - x[0]) / 50000.0) / 31536000.0"), V_cg)
File(out_dir + "u_b.xml") << u_b
File(out_dir + "u_b.pvd") << u_b


### Create bed and surface functions

# Maximum ice thickness
h_max = 1500.
# Length of ice sheet 
length = 50e3
# Center of trough 
center = 10e3
# Maximum trough depth
depth = 200.0
      
class Bed(Expression):
  def eval(self,value,x):
    value[0] = sin(pi / 180) * x[0]

class Surface(Expression):
  def eval(self,value,x):
    value[0] = sqrt((x[0] + 50.0) * h_max**2 / length)

# Surface
S = project(Surface(), V_cg)
# Bed elevation
B = project(Bed(), V_cg)
# Ice thickness
H = project(S - B, V_cg)

File(out_dir + "B.xml") << B
File(out_dir + "B.pvd") << B

File(out_dir + "H.xml") << H
File(out_dir + "H.pvd") << H

rho_i = pcs['rho_i']
rho_w = pcs['rho_w']
g = pcs['g']

# Potential at 0 pressure
phi_m = project(rho_w * g * B, V_cg)
# Overburden pressure
p_i = project(rho_i * g * H, V_cg)
# Potential at overburden pressure
phi_0 = project(phi_m + p_i, V_cg)

File(out_dir + "phi_m.xml") << phi_m
File(out_dir + "phi_m.pvd") << phi_m
File(out_dir + "p_i.xml") << p_i
File(out_dir + "p_i.pvd") << p_i
File(out_dir + "phi_0.xml") << phi_0
File(out_dir + "phi_0.pvd") << phi_0


### Create a facet function with marked boundaries

# Margin
class MarginSub(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], 0.0)
    
# Divide 
class DivideSub(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], 50000.0)
    
# Near margin channels
class NearMargin(SubDomain):
  def inside(self, x, on_boundary):
    return not on_boundary and x[0] < 1e3 and abs(x[1] - 10e3) < 3500.0
    
nm = NearMargin()  
ms = MarginSub()
ds = DivideSub()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
ms.mark(boundaries, 1)
ds.mark(boundaries, 2)
nm.mark(boundaries, 3)

File(out_dir + "boundaries.xml") << boundaries
File(out_dir + "boundaries.pvd") << boundaries
plot(boundaries, interactive = True)


### Create a mask for the ODE solver

# Create a mask array that is 0 on mesh edges and 1 on interior edges. This 
# is used to prevent opening on exterior edges
v_cr = TestFunction(V_cr)
mask = assemble(v_cr('+') * dS)
mask[mask.array() > 0.0] = 1.0

mask_cr = Function(V_cr)
mask_cr.vector()[:] = mask.array()

File(out_dir + "mask.xml") << mask_cr
File(out_dir + "mask.pvd") << mask_cr