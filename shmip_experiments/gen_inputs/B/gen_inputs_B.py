"""
Model inputs for intercomparisson tests B1 - B6. 
"""

from dolfin import *
from constants import *
import numpy as np

ns = range(1,6)

out_files = ["../../inputs/B/input_B" + str(n) + ".hdf5" for n in ns]


# Directory to write model inputs
mesh = Mesh("../../inputs/mesh/mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

# Sliding speed
u_b = interpolate(Constant(1e-6), V_cg)
# Conductivity 
k = interpolate(Constant(5e-3), V_cg)
# Initial sheet height  
h = interpolate(Constant(0.09), V_cg)

# Length of ice sheet 
length = 100e3
      
class Bed(Expression):
  def eval(self,value,x):
    value[0] = 0.0

class Surface(Expression):
  def eval(self,value,x):
    value[0] = 6.0 * (sqrt(x[0] + 5e3) - sqrt(5e3)) + 1

# Surface
S = project(Surface(degree=1), V_cg)
# Bed elevation
B = project(Bed(degree=1), V_cg)
# Ice thickness
H = project(S - B, V_cg)

# Potential at 0 pressure
phi_m = project(pcs['rho_w'] * pcs['g'] * B, V_cg)
# Ice overburden pressure
p_i = project(pcs['rho_i'] * pcs['g'] * H, V_cg)
# Potential at overburden pressure
phi_0 = project(phi_m + p_i, V_cg)


### Create a facet function with marked boundaries

# Margin
class MarginSub(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], 0.0)
    
ms = MarginSub()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
ms.mark(boundaries, 1)


# Get x and y coordinates of mesh vertices
coords = V_cg.tabulate_dof_coordinates().reshape(V_cg.dim(), 2)
coords_x = coords[:,0]
coords_y = coords[:,1]



for n in range(len(ns)):
  
  ## Calculate the melt function 
  
  # Moulin locations and fluxes
  data = np.loadtxt('B' + str(ns[n]) + '_M.csv', delimiter=',')
  # Moulin melt function
  m = Function(V_cg)
  # Initial condition for h
  h = Function(V_cg)
  
  if n == 0:
      data = np.array([data])
  
  # Moulin x locations
  moulin_xs = data[:,1]
  # Moulin y locations
  moulin_ys = data[:,2]
  # Flux into moulin in m^3 / s
  fluxes = data[:,3]
  # DOF indices of moulin vertexes
  indexes = np.zeros(len(moulin_xs))
  # DOF values at moulin vertexes
  values = np.zeros(len(moulin_xs))
  
  
  # Find moulin dofs and values
  for i in range(len(moulin_xs)):  
    moulin_x = moulin_xs[i]
    moulin_y = moulin_ys[i]
    flux = fluxes[i]
    
    # Find the vertex nearest to each moulin
    distances = np.sqrt((coords_x - moulin_x)**2 + (coords_y - moulin_y)**2)
    dof_index = distances.argmin()
    
    # Make a point source function and integrate it
    m.vector()[dof_index] = 1
    int_m = assemble(m * dx)
    indexes[i] = dof_index
    values[i] = flux / int_m
  
    # Reset m to zero
    m.vector()[dof_index] = 0.0
    
  # Create the melt function
  m.vector()[:] = 7.93e-11
  m.vector()[indexes] += values
  h.vector()[:] = 0.09
  h.vector()[indexes] = 1.0
  
  
  ## Write input file 
  
  # Write inputs to a hdf5 file
  f = HDF5File(mesh.mpi_comm(), out_files[n], 'w')
  f.write(mesh, "mesh")
  f.write(u_b, "u_b_0")
  f.write(k, "k_0")
  f.write(B, "B")
  f.write(H, "H")
  f.write(phi_0, "phi_0")
  f.write(boundaries, "boundaries")  
  f.write(m, "m_0")
  f.write(h, "h_0")
  
  f.close()
