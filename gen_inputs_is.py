"""
This script generates all fields such as ice thickness, bed elevation, and melt
needed to run the model for an ice sheet with a sloped bed.
"""

from dolfin import *
from constants import *

# Input directory
in_dir = "inputs_is/"
# Directory to write model inputs
out_dir = "inputs_is/"

mesh = Mesh("inputs_is/mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Bed
B = Function(V_cg)
File(in_dir + "B.xml") >> B

# Thickness
H = Function(V_cg)
File(in_dir + "H.xml") >> H

# Constants
rho_i = pcs['rho_i']
rho_w = pcs['rho_w']
g = pcs['g']

# Potential at 0 pressure
phi_m = project(rho_w * g * B, V_cg)
# Overburden pressure
p_i = project(rho_i * g * H, V_cg)
# Potential at overburden
phi_0 = project(phi_m + p_i, V_cg)

File(out_dir + "phi_m.xml") << phi_m
File(out_dir + "phi_m.pvd") << phi_m
File(out_dir + "p_i.xml") << p_i
File(out_dir + "p_i.pvd") << p_i
File(out_dir + "phi_0.xml") << phi_0
File(out_dir + "phi_0.pvd") << phi_0
