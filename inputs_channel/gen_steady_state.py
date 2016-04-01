

from dolfin import *


input_file = HDF5File(mpi_comm_world(), 'out.hdf5', 'r')
output_file = HDF5File(mpi_comm_world(), 'channel_ref_steady.hdf5', 'w')

# Load the mesh
mesh = Mesh()
input_file.read(mesh, "mesh", False)    
V_cg = FunctionSpace(mesh, "CG", 1) 
V_cr = FunctionSpace(mesh, "CR", 1) 


### Read variables

# Load the most recent cavity height value
h = Function(V_cg)
num_steps = input_file.attributes("h")['count']
h_last = "h/vector_" + str(num_steps - 1)
input_file.read(h, h_last)

# Most recent channel height
S = Function(V_cr)
S_last = "S/vector_" + str(num_steps - 1)
input_file.read(S, S_last)

# Most recent potential
phi = Function(V_cr)
phi_last = "phi/vector_" + str(num_steps - 1)
input_file.read(phi, phi_last)


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

# Edge lengths
edge_lens = Function(V_cr)
input_file.read(edge_lens, 'edge_lens')

# Sheet conductivity
k_0 = Function(V_cg)
input_file.read(k_0, 'k_0')

# Channel conductivity
k_c_0 = Function(V_cg)
input_file.read(k_c_0, 'k_c_0')

# Mask
mask = Function(V_cr)
input_file.read(mask, 'mask')


### Write variables

output_file.write(mesh, "mesh")
output_file.write(B, "B")
output_file.write(H, "H")
output_file.write(m, 'm_0')
output_file.write(u_b, 'u_b_0')
output_file.write(h, "h_0")
output_file.write(S, "S_0")
output_file.write(boundaries, "boundaries")
output_file.write(edge_lens, "edge_lens")
output_file.write(k_0, "k_0")
output_file.write(k_c_0, "k_c_0")
output_file.write(mask, "mask")
output_file.write(phi, "phi_0")