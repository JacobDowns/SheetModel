"""
This script generates all fields such as ice thickness, bed elevation, and melt
needed to run the model for an ice sheet with a sloped bed.
"""

from dolfin import *
from constants import *
import numpy as np

# Directory to write model inputs
out_dir = "inputs_trough/"

mesh = Mesh(out_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Melt
m = project(Expression("(1.0 + 1.5 * (50000.0 - x[0]) / 50000.0) / 31536000.0"), V_cg)
File(out_dir + "m.xml") << m
File(out_dir + "m.pvd") << m

# Sliding speed
u_b = project(Expression("(10.0 + 240.0 * (50000.0 - x[0]) / 50000.0) / 31536000.0"), V_cg)
File(out_dir + "u_b.xml") << u_b
File(out_dir + "u_b.pvd") << u_b

# Bed and surface functions

# Maximum ice thickness
h_max = 1500.
# Length of ice sheet 
length = 50e3
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
      value[0] = -depth * e**(-((Y - center) / (2500.0 * shape(X, 500.0)) )**2.0) * shape(X, 7000.0)

class Surface(Expression):
  def eval(self,value,x):
    value[0] = sqrt((x[0] + 20.0) * h_max**2 / length)

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
    
ms = MarginSub()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
ms.mark(boundaries, 1)

File(out_dir + "boundaries.xml") << boundaries
File(out_dir + "boundaries.pvd") << boundaries
plot(boundaries, interactive = True)

