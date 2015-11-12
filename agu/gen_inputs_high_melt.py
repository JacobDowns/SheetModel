"""
This script generates all fields such as ice thickness, bed elevation, and melt
needed to run the model for an ice sheet on a trough. This script generates
a high melt rate. 
"""

import sys
from dolfin import *
import numpy as np
sys.path.insert(0, '/home/jake/schoof/')
from constants import *

# Directory to write model inputs
out_dir = "inputs_high_melt/"
mesh = Mesh(out_dir + "mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Melt
m = project(Expression("(1.0 + (4.0 * (60000.0 - x[0]) / 60000.0)) / 31536000.0"), V_cg)
File(out_dir + "m.xml") << m
File(out_dir + "m.pvd") << m

# Sliding speed
u_b = project(Expression("(50.0 + 250.0 * (60000.0 - x[0]) / 60000.0) / 31536000.0"), V_cg)
File(out_dir + "u_b.xml") << u_b
File(out_dir + "u_b.pvd") << u_b

# Bed and surface functions

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

