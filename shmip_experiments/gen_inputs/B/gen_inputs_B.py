"""
Model inputs for intercomparisson tests B1 - B6. 
"""

from dolfin import *
from constants import *
import numpy as np



# Directory to write model inputs
mesh = Mesh("../../inputs/mesh/mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

#ns = range(1,6)
ns = [1]

for n in ns:
  data = np.loadtxt('B' + str(n) + '_M.csv', delimiter=',')
  m = Function(V_cg)
  
  moulin_xs = data[:,1]
  moulin_ys = data[:,2]
  flux = data[:,3]
  
  for i in range(moulin_xs):
    moulin_x = moulin_xs[i]
    moulin_y = moulin_ys[i]
    
    distances = np.sqrt((coords_x - xd)**2 + (coords_y - yd)**2)
    i = distances.argmin()
    
    
  
  
quit()
                              



out_files = ["../../inputs/A/input_A" + str(n) + ".hdf5" for n in ns]
melt_rates = [7.93e-11, 1.59e-9, 5.79e-9, 2.5e-8, 4.5e-8, 5.79e-7]
sheet_heights = [0.005, 0.001, 0.0015, 0.002, 0.0025, 0.003]



coords = V_cg.tabulate_dof_coordinates().reshape(V_cg.dim(), 2)
coords_x = coords[:,0]
coords_y = coords[:,1]

xd = 59000.0
yd = 8000.0 




f = Function(V_cg)
f.vector()[i] = 1.0

File('f.pvd') << f


quit()

# Sliding speed
u_b = interpolate(Constant(1e-6), V_cg)
# Conductivity 
k = interpolate(Constant(5e-3), V_cg)

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

# Initial sheet height


for n in range(6):
  # Write inputs to a hdf5 file
  f = HDF5File(mesh.mpi_comm(), out_files[n], 'w')
  f.write(mesh, "mesh")
  f.write(u_b, "u_b_0")
  f.write(k, "k_0")
  f.write(B, "B")
  f.write(H, "H")
  f.write(phi_0, "phi_0")
  f.write(boundaries, "boundaries")  

  # Get initial sheet height  
  h = interpolate(Constant(sheet_heights[n]), V_cg)
  f.write(h, "h_0")

  # Get melt rate
  m = interpolate(Constant(melt_rates[n]), V_cg)
  f.write(m, "m_0")
  
  f.close()
  
 

