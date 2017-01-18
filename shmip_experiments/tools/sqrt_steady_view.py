# -*- coding: utf-8 -*-
"""
Loads an hdf5 results file and writes output as netcdf.
"""

from steady_view import *

class SqrtSteadyView(SteadyView):
  
  def __init__(self, input_file):
    SteadyView.__init__(self, input_file)
    
    ## Setup some stuff for doing width and length-wise integration across the domain    
    
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
    self.length_bc = DirichletBC(self.V_cg, 0.0, boundaries, 2)
    self.width_bc = DirichletBC(self.V_cg, 0.0, boundaries, 1)
    
    # Test function
    self.theta = TestFunction(self.V_cg)
    
 
  # Width integrate a field 
  def width_integrate(self, u):
    v = TrialFunction(self.V_cg)
    a = v.dx(1) * self.theta * dx
    L = u * self.theta * dx
    v = Function(self.V_cg)
    solve(a == L, v, self.width_bc)
    return v
    
    
  # Width integrate effective pressure
  def width_integrate_N(self):
    return self.width_integrate(self.N)
    
    
  # Width integrate sheet height
  def width_integrate_h(self):
    return self.width_integrate(self.h)
    
  
  # Width integrate q
  def width_integrate_q(self):
    return self.width_integrate(self.q)
    