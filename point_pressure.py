# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:00:39 2015

@author: jake
"""

from dolfin import *

# Model input directory
in_dir = "inputs_trough_60_20/"
# Output directory
out_dir = "out_trough_time_60_20/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"

# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

pfo = Function(V_cg)
for i in range(548):
  File(check_dir + "pfo_" + str(i) + ".xml") >> pfo
  
  # Sample the pressure at a few places
  p_10 = pfo([0.0])