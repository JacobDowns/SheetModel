"""
Load all of the Issunguata model inputs and stick them in an hdf5 file. Sheet
thickness is the steady state sheet thickness.
"""

from dolfin import *

# Directory to write model inputs
out_dir = "inputs_is/"
mesh = Mesh(out_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Write inputs to a hdf5 file
f = HDF5File(mesh.mpi_comm(), out_dir + "inputs_is_steady.hdf5", 'w')

# Write the mesh to a file
f.write(mesh, "mesh")

# Bed
B = Function(V_cg)
File(out_dir + "B.xml") >> B
f.write(B, "B")

# Ice thickness
H = Function(V_cg)
File(out_dir + "H.xml") >> H
f.write(H, "H")

# Melt
m = Function(V_cg)
File(out_dir + "m.xml") >> m
f.write(m, "m_0")

# Sliding speed
u_b = Function(V_cg)
File(out_dir + "u_b.xml") >> u_b
f.write(u_b, "u_b_0")

# Facet function marking boundaries
boundaries = FacetFunction("size_t", mesh)
File(out_dir + "boundaries.xml") >> boundaries
f.write(boundaries, "boundaries")

# Initial sheet height
h = Function(V_cg)
File(out_dir + "h_99.xml") >> h
f.write(h, "h_0")