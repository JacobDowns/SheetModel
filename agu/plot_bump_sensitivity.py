# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:45:41 2015

@author: jake
"""

from dolfin import *
from pylab import *
#from constants import *

# Model input directory
in_dir = "inputs_flat_bed/"

# Process number
MPI_rank = MPI.rank(mpi_comm_world())

# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Pressures for each smiulation
pfo = Function(V_cg)

# Points to record pressure
xs = [10e3, 20e3, 50e3]
ys = [10e3, 10e3, 10e3]

# List that contains lists of pressures for each point 
pfos = [[] for i in range(len(xs))]


### Load all the data
for i in range(1, 11):
  print ("i", i)
  
  # Load overburden pressure
  File("out_bump_sensitivity/pfo" + str(i) + ".xml") >> pfo
  
  for i in range(len(xs)):
    x = xs[i]
    y = ys[i]
    
    # Record pressures
    pfos[i].append(pfo([x, y]))
    
print pfos[0]
    
h_rs = arange(0.06, 0.16, 0.01)[::-1]
print h_rs
    
plot(h_rs, pfos[0], 'ro-')
plot(h_rs, pfos[1], 'go-')
plot(h_rs, pfos[2], 'bo-')
show()

 