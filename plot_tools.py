# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 13:08:33 2016

@author: jake
"""

from dolfin import *
from constants import *
from cr_tools import *

class PlotTools(object):
  
  def __init__(self, input_file):

    # Load a results file
    self.input_file = HDF5File(mpi_comm_world(), input_file, 'r') 

    # Load the mesh
    self.mesh = Mesh()
    self.input_file.read(self.mesh, "mesh", False)  
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)
    self.V_cr = FunctionSpace(self.mesh, "CR", 1)
    
    # Create a CR tools object for plotting CR functions
    self.edge_lens = Function(self.V_cr)
    self.input_file.read(self.edge_lens, "edge_lens")
    self.cr_tools = CRTools(self.mesh, self.V_cg, self.V_cr, self.edge_lens)
    
    # Get the number of time steps
    self.num_steps = 0
    try:
      self.num_steps = self.input_file.attributes("h")['count']
    except:
      pass
    
    # Compute overburden potential
    B = Function(self.V_cg)
    H = Function(self.V_cg)
    self.input_file.read(B, "B")
    self.input_file.read(H, "H")
    # Potential at 0 pressure
    phi_m = project(pcs['rho_w'] * pcs['g'] * B, self.V_cg)
    # Ice overburden pressure
    p_i = project(pcs['rho_i'] * pcs['g'] * H, self.V_cg)
    # Potential at overburden pressure
    self.phi0 = project(phi_m + p_i, self.V_cg)
    
    # Compute the total area of the domain
    f = interpolate(Constant(1.0), self.V_cg)
    self.area = assemble(f * dx)
    
    # Initial conductivity
    self.k_0 = Function(self.V_cg)
    self.input_file.read(self.k_0, "k_0")
    # Initial melt
    self.m_0 = Function(self.V_cg)
    self.input_file.read(self.m_0, "m_0")    
    # Initial potential 
    self.phi_0 = Function(self.V_cg)
    self.input_file.read(self.phi_0, "phi_0")
    
    # Fields we might wish to plot
    
    # Channel area
    self.S = Function(self.V_cr)
    #Channel area as a facet function
    self.S_ff = FacetFunctionDouble(self.mesh)
    # Sheet height
    self.h = Function(self.V_cg)
    # Hydraulic potential
    self.phi = Function(self.V_cg)
    # Presure as a fraction of overburden
    self.pfo = Function(self.V_cg)
    # Melt 
    self.m = Function(self.V_cg)
    # Conductivity
    self.k = Function(self.V_cg)
    

  # Get the initial potential 
  def get_phi_0(self):
    return self.phi_0


  # Get S at the i-th time step
  def get_phi(self, i):
    if i < self.num_steps:
      self.input_file.read(self.phi, "phi/vector_" + str(i))
      return self.phi
        
  
  # Get h at the i-th time step
  def get_h(self, i):
    if i < self.num_steps:
      self.input_file.read(self.h, "h/vector_" + str(i))
      return self.h
        
        
  # Get pfo at the i-th time step
  def get_pfo(self, i):
    if i < self.num_steps:
      self.input_file.read(self.pfo, "phi/vector_" + str(i))
      self.pfo.vector()[:] = self.pfo.vector().array() / self.phi0.vector().array()
      return self.pfo
        
        
  # Compute the average pfo over the domain
  def get_avg_pfo(self, i):
    if i < self.num_steps:
      self.get_pfo(i)
      return assemble(self.pfo * dx) / self.area
     
    
  # Get S at the i-th time step
  def get_S(self, i):
    if i < self.num_steps:
      self.input_file.read(self.S, "S/vector_" + str(i))
      return self.S


  # Get S at the i-th time step as a facet function
  def get_S_ff(self, i):
    if i < self.num_steps:
      self.get_S(i)
      self.cr_tools.copy_cr_to_facet(self.S, self.S_ff)
      return self.S_ff
        
  
  # Get m at the i-th time step
  def get_m(self, i):
    if i < self.num_steps:
      try:
        self.input_file.read(self.m, "m/vector_" + str(i))
        return self.m
      except:
        self.m.assign(m_0)
        return self.m
        
        
  # Get the total melt input at the i-th time step
  def get_total_m(self, i):
    if i < self.num_steps:
      self.get_m(i)
      return assemble(self.m * dx)
        
  
  # Get k at the i-th time step
  def get_k(self, i):
    if i < self.num_steps:
      try:
        self.input_file.read(self.k, "k/vector_" + str(i))
        return self.k
      except:
        self.k.assign(self.k_0)
        return self.k
        
        
  # Compute the volume of water in the sheet at the i-th time
  def get_sheet_volume(self, i):
    if i < self.num_steps:
      self.get_h(i)
      return assemble(self.h * dx)
      
  
  # Writes h, S, and pfo out to pvd files for each time step
  def write_pvds(self, out_file):
    out_h = File(out_file + 'h.pvd')
    out_pfo = File(out_file + 'pfo.pvd')
    out_S = File(out_file + 'S.pvd')
    
    for i in range(self.num_steps):
      print ("i", i)
      out_h << self.get_h(i)
      out_pfo << self.get_pfo(i)
      out_S << self.get_S_ff(i)
      
      