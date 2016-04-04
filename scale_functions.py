# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:42:36 2016
@author: jake
"""
from dolfin import *
from constants import *

class ScaleFunctions(object):
  
  def __init__(self, input_file, k_min = 5e-5, k_max = pcs['k'], m_max = 5.0, shutoff_length = 30.0 * pcs['spd'], u_b_max = 100.0, lag_time = 0.0):
    f = HDF5File(mpi_comm_world(), input_file, 'r') 
    
    # Get the mesh
    mesh = Mesh()
    f.read(mesh, "mesh", False)  
    self.V_cg = FunctionSpace(mesh, "CG", 1)

    # Get melt and sliding speed
    self.m = Function(self.V_cg)
    self.u_b = Function(self.V_cg)
    f.read(self.m, "m_0")
    f.read(self.u_b, "u_b_0")
       
    # Shutoff length
    self.shutoff_length = shutoff_length
    # Minimum conductivity
    self.k_min = k_min
    # Maximum conductivity
    self.k_max = k_max
    # Scaling parameter that sets the maximum possible conductivity
    self.a = (k_max - k_min) / (m_max / pcs['spy'])
    # Lag of conductivity behind melt
    self.b = lag_time
    # Parameter in the sliding speed scale function 
    self.c = (u_b_max / pcs['spy']) / self.u_b.vector().max()

    
  # Melt scale function
  def m_scale(self, t):
    if t < 0.0:
      return 1.0
    if t <= self.shutoff_length:
      return cos((pi / (2.0 * self.shutoff_length)) * t)
    return 0.0
  
  
  # Conductivity scale function
  def k_scale(self, t):
    return self.a * self.m_scale(t - self.b)


  # Sliding speed scale function
  def u_b_scale(self, t):
    return -self.m_scale(t) * (self.c - 1.0) + self.c

    
  # Gets the melt function at a particular time
  def get_m(self, t):
    return project(Constant(self.m_scale(t)) * self.m, self.V_cg)


  # Gets the conductivity function at a particular time
  def get_k(self, t):
    return project(Constant(self.k_scale(t)) * self.m + Constant(self.k_min), self.V_cg)
    
    
  # Gets the sliding speed at a particular time
  def get_u_b(self, t):
    return project(Constant(self.u_b_scale(t)) * self.u_b, self.V_cg)
    