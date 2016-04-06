from dolfin import *
from constants import *
from plot_tools import *
import numpy as np

""" An oddly specific class specifically for visualizing various channel model
fields on a rectangular domain. """

class PlotRect(PlotTools):
  
  def __init__(self, input_file):
    PlotTools.__init__(self, input_file)
    
    ### Setup some stuff for doing width and length-wise integration across the domain    
    
    # Mark the bottom edge and right edge of the mesh with a facet function
    class BottomEdge(SubDomain):
      def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0)
        
    class LeftEdge(SubDomain):
      def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)
        
    be = BottomEdge()
    le = LeftEdge()
    boundaries = FacetFunction("size_t", self.mesh)
    boundaries.set_all(0)
    be.mark(boundaries, 1)
    le.mark(boundaries, 2)
    
    # Boundaries for doing length and width integration
    self.x_bc = DirichletBC(self.V_cg, 0.0, boundaries, 2)
    self.y_bc = DirichletBC(self.V_cg, 0.0, boundaries, 1)
    
    # Test function
    self.w = TestFunction(self.V_cg)
    
    
    ### Expression for sheet discharge
    
    phi_reg = Constant(1e-15)
    alpha = pcs['alpha']
    delta = pcs['delta']
    
    q = -self.k * self.h**alpha * (dot(grad(self.phi), grad(self.phi)) + phi_reg)**(delta / 2.0) * grad(self.phi)
    self.qn = dot(q, as_vector([-1, 0]))
    
    increments = 250
    self.xs = np.linspace(1, 60e3, increments)
    
    
    ### Compute some stuff necessary for computing the integrated channel discharge
    
    self.edge_to_dof1, self.edge_to_dof2 = self.calculate_edge_to_dof_maps()

    # Get the x coordinates of the start end end points of each edge
    xs_func = project(Expression("x[0]"), self.V_cg)
    edge_x1 = Function(self.V_cr)
    edge_x2 = Function(self.V_cr)
    edge_x1.vector()[:] = xs_func.vector().array()[self.edge_to_dof1] 
    edge_x2.vector()[:] = xs_func.vector().array()[self.edge_to_dof2] 
    self.edge_x1 = edge_x1.vector().array()
    self.edge_x2 = edge_x2.vector().array()
    
    # Derivative of potential along channel edges
    self.dphi_ds = Function(self.V_cr)
    # Channel conductivity
    k_c = Function(self.V_cg)
    self.input_file.read(k_c, "k_c_0")
    self.k_c = k_c.vector().array()[0]
    self.delta = delta
    self.alpha = alpha
    
    
  # Integrate a field across the width of the ice sheet
  def integrate_y(self, u):
    v = TrialFunction(self.V_cg)
    a = v.dx(1) * self.w * dx
    L = u * self.w * dx
    v = Function(self.V_cg)
    solve(a == L, v, self.y_bc)
    return v
    
    
  # Integrate a field across the length of the ice sheet
  def integrate_x(self, u):
    v = TrialFunction(self.V_cg)
    a = v.dx(0) * self.w * dx
    L = u * self.w * dx
    v = Function(self.V_cg)
    solve(a == L, v, self.x_bc)
    return v
    
  
  # Get the width integrated sheet discharge at the ith time step
  def get_int_sheet_discharge(self, i):
    if i < self.num_steps:
      self.get_h(i)
      self.get_phi(i)
      self.get_k(i)
      
      int_qn = self.integrate_y(self.qn)
      
      qs = []
      for x in self.xs:
        qs.append(int_qn([x, 20e3]))
        
      return np.array(qs)
      
  # Get the width integrated sheet discharge at the ith time step
  def get_int_sheet_height(self, i):
    if i < self.num_steps:
      self.get_h(i)
      
      int_h = self.integrate_y(self.h)
      
      hs = []
      for x in self.xs:
        hs.append(int_h([x, 20e3]))
        
      return np.array(hs) / 20e3
      
    
  # Get the width integrated channel discharge at the ith time step across a line at x
  def get_channel_discharge_x(self, i, x):  
    if i < self.num_steps:
      self.get_phi(i)
      self.get_S(i) 
  
      #phi_s = self.compute_dphi_ds(i)
      self.cr_tools.ds(self.phi, self.dphi_ds)
      phi_s = self.dphi_ds.vector().array()     
      S_n = self.S.vector().array()
      
      phi_reg = 1e-15
      Q_n = self.k_c * S_n**self.alpha * abs(phi_s + phi_reg)**self.delta * np.absolute(phi_s)
      
      c1 = np.logical_and(self.edge_x1 <= x, self.edge_x2 >= x)
      c2 = np.logical_and(self.edge_x1 >= x, self.edge_x2 <= x)
      indexes = np.logical_or(c1, c2)
    
      discharge = sum(Q_n[indexes])
      return np.array(discharge)
      
      
  # Compute oriented edge derivatives      
  def compute_dphi_ds(self, i):
    self.edge_x1
    self.edge_x2
    
    phi_v1 = self.phi.vector().array()[self.edge_to_dof1] 
    phi_v2 = self.phi.vector().array()[self.edge_to_dof2] 
    
    phi_s = np.zeros(len(self.edge_x1))  
    
    indexes = self.edge_x1 > self.edge_x2

    phi_s[indexes] = phi_v1[indexes] - phi_v2[indexes]
    indexes = np.logical_not(indexes)
    phi_s[indexes] = phi_v2[indexes] - phi_v1[indexes]
    
    phi_s /= self.edge_lens.vector().array()
    
    return phi_s
    
        
  # Get the width integrated channel discharge at the ith time step
  def get_int_channel_discharge(self, i):  
    if i < self.num_steps:
      Qs = []
      for x in self.xs:
        Qs.append(self.get_channel_discharge_x(i,x))
        
      return Qs
    
        
  # This calculates the mapping from facet dof indices to facets.  It is
  # analogous to the V.dofmap().dof_to_vertex_map(mesh) method.
  def calculate_edge_to_facet_map(self):
    n_V = self.V_cr.dim()
  
    # Find coordinates of dofs and put into array with index
    coords_V = np.hstack((np.reshape(self.V_cr.dofmap().tabulate_all_coordinates(self.mesh),(n_V,2)), np.zeros((n_V,1))))
    coords_V[:,2] = range(n_V)
  
    # Find coordinates of facets and put into array with index
    coords_f = np.zeros((n_V,3))
    for f in dolfin.facets(self.mesh):
        coords_f[f.index(),0] = f.midpoint().x()
        coords_f[f.index(),1] = f.midpoint().y()
        coords_f[f.index(),2] = f.index() 
  
    # Sort these the same way
    coords_V = np.array(sorted(coords_V,key=tuple))
    coords_f = np.array(sorted(coords_f,key=tuple))
  
    # the order of the indices becomes the map
    V2fmapping = np.zeros((n_V,2))
    V2fmapping[:,0] = coords_V[:,2]
    V2fmapping[:,1] = coords_f[:,2]
  
    return (V2fmapping[V2fmapping[:,0].argsort()][:,1]).astype('int')
  
  
  # Computes maps from each edge in a CR function space to the associated
  # vertex dofs on the edge    
  def calculate_edge_to_dof_maps(self):
    # First vertex index associated with each facet
    f_0 = dolfin.FacetFunction('uint', self.mesh)    
    # Second vertex index associated with each facet
    f_1 = dolfin.FacetFunction('uint', self.mesh)
  
    # Map from vertex index to degree of freedom
    v2d = vertex_to_dof_map(self.V_cg)
    
    # Get the two dof indexes associated with each facet
    for f in dolfin.facets(self.mesh):
      # Vertexes associated with this facet
      v0, v1 = f.entities(0)
      
      # The first dof associated with this facet
      f_0[f] = v2d[v0]
      # The second dof associated with this facet
      f_1[f] = v2d[v1]
    
    edge_to_facet_map = self.calculate_edge_to_facet_map()
    edge_to_dof0 = f_0.array()[edge_to_facet_map]
    edge_to_dof1 = f_1.array()[edge_to_facet_map]
    
    return (edge_to_dof0, edge_to_dof1)
  

        
      
      