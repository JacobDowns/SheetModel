from dolfin import *
import numpy as np

# Computes directional of CG functions along edges as well as midpoints of CG
# functions along edges in parallel
class CRTools(object):
  
  def __init__(self, mesh, V_cg, V_cr, edge_lens) :
    self.mesh = mesh
    # DG function Space
    self.V_cg = V_cg
    # CR function space
    self.V_cr = V_cr
    # Process
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # Compute a map from local facets to global edge indexs
    self.compute_local_facet_to_global_edge_index_map()
    # Compute a map from local edges to global edge indexes
    self.compute_local_edge_to_global_edge_index_map()
    # Compute local edge lengths
    self.edge_lens = edge_lens
    
    # We'll set up a form that allows us to take the derivative of CG functions
    # over edges 
    self.U = Function(V_cg)
    # CR test function
    v_cr = TestFunction(V_cr)
    # Facet and tangent normals
    n = FacetNormal(mesh)
    t = as_vector([n[1], -n[0]])
    # Directional derivative form
    self.F = (dot(grad(self.U), t) * v_cr)('+') * dS
    
    # Facet function for plotting 
    self.ff_plot = FacetFunctionDouble(mesh)
    
    
  # Copies a CR function to a facet function
  def copy_cr_to_facet(self, cr, ff) :
    # Gather all edge values from each of the local arrays on each process
    cr_vals = Vector()
    cr.vector().gather(cr_vals, np.array(range(self.V_cr.dim()), dtype = 'intc'))
    # Get the edge values corresponding to each local facet
    local_vals = cr_vals[self.local_facet_to_global_edge_index_map]    
    ff.array()[:] = local_vals

  
  # Compute a map from local facets to indexes in a global cr vector
  def compute_local_facet_to_global_edge_index_map(self):
    # Compute the midpoints of local edges in a cr function
    edge_coords = self.V_cr.dofmap().tabulate_all_coordinates(self.mesh)
    edge_coords = edge_coords.reshape((len(edge_coords)/2), 2)
    
    # Gather a global array of midpoints for a cr function
    fx = Function(self.V_cr)
    fy = Function(self.V_cr)
    fx.vector().set_local(edge_coords[:,0])
    fy.vector().set_local(edge_coords[:,1])
    fx.vector().apply("insert")
    fy.vector().apply("insert")
    vecx = Vector()
    vecy = Vector()    
    fx.vector().gather(vecx, np.array(range(self.V_cr.dim()), dtype = 'intc'))
    fy.vector().gather(vecy, np.array(range(self.V_cr.dim()), dtype = 'intc'))
    
    # x coordinate of midpoint
    global_edge_coords_x = vecx.array()
    # y coordinate of midpoint
    global_edge_coords_y = vecy.array()
    
    # Create a dictionary that maps a midpoint tuple to an index in the global cr vector
    midpoint_to_cr_index = {}
    
    for i in range(len(global_edge_coords_x)):
      x = global_edge_coords_x[i]
      y = global_edge_coords_y[i]
      midpoint_to_cr_index[(x,y)] = i
      
    # Create an array that maps local facet indexes to indexes in the global cr vector
    local_facet_to_global_cr = []
    i = 0
    for f in dolfin.facets(self.mesh):
      x = f.midpoint().x()
      y = f.midpoint().y()
      v0, v1 = f.entities(0)
      
      local_facet_to_global_cr.append(midpoint_to_cr_index[(x,y)])
      i += 1
      
    self.local_facet_to_global_edge_index_map = np.array(local_facet_to_global_cr)
    
    
  # Computes a map from local edges to to indexes in a global cr vector
  def compute_local_edge_to_global_edge_index_map(self):
    f = Function(self.V_cr)
    self.local_edge_to_global_edge_index_map = np.zeros(len(f.vector().array()), dtype = 'intc')
    for i in range(len(f.vector().array())):
      self.local_edge_to_global_edge_index_map[i] = self.V_cr.dofmap().local_to_global_index(i)

  # Computes the directional derivative of a CG function over edges 
  def ds(self, cg, cr):
    self.U.assign(cg)
    
    # Get the height difference of two vertexes on each edge
    A = abs(assemble(self.F).array())
    # Now divide by the edge lens
    dcg_ds = A / self.edge_lens.vector().array()
    
    cr.vector().set_local(dcg_ds)
    cr.vector().apply("insert") 

  
  # Computes the value of a CG functions at the midpoint of edges and copies
  # the result to a CR function
  def midpoint(self, cg, cr):
    cr.assign(interpolate(cg, self.V_cr))
  
  
  # Plots a CR function
  def plot_cr(self, cr):
    self.copy_cr_to_facet(cr, self.ff_plot)
    plot(self.ff_plot, interactive = True)
