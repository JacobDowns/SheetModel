# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:59:04 2016

@author: jake
"""

from dolfin import *


input_file = HDF5File(mpi_comm_world(), 'inputs_is_outlets_steady.hdf5', 'r')
output_file = HDF5File(mpi_comm_world(), 'inputs_is_outlets_steady1.hdf5', 'w')

# Load the mesh
mesh = Mesh()
input_file.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1) 


### Read variables

# Load the most recent cavity height value
h = Function(V_cg)
num_steps = input_file.attributes("h")['count']
h_last = "h/vector_" + str(num_steps - 1)
input_file.read(h, h_last)

# Bed
B = Function(V_cg)
input_file.read(B, 'B')

# Ice thickness
H = Function(V_cg)
input_file.read(H, 'H')

# Melt
m = Function(V_cg)
input_file.read(m, 'm_0')

# Sliding
u_b = Function(V_cg)
input_file.read(u_b, 'u_b_0')

# Boundaries
boundaries = FacetFunction('size_t', mesh)
input_file.read(boundaries, 'boundaries')


### Write variables

output_file.write(mesh, "mesh")
output_file.write(B, "B")
output_file.write(H, "H")
output_file.write(m, 'm_0')
output_file.write(u_b, 'u_b_0')
output_file.write(h, "h_0")
output_file.write(boundaries, "boundaries")
