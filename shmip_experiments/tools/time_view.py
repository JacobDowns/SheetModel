# -*- coding: utf-8 -*-
"""
Loads an hdf5 results file. Writes a netcdf file. 
"""

from dolfin import *
from constants import *
from netCDF4 import Dataset
import numpy as np

class TimeView(object):
  
  def __init__(self, input_file):

    # Load a results file
    self.input_file = HDF5File(mpi_comm_world(), input_file, 'r') 

    # Load the mesh
    self.mesh = Mesh()
    self.input_file.read(self.mesh, "mesh", False)  
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)
    
    
    # Get the number of time steps
    self.num_steps = 0
    try:
      self.num_steps = self.input_file.attributes("h")['count']
    except:
      raise Exception("Not a results file.")
      sys.exit(1)
    
    # Bed
    self.B = Function(self.V_cg)
    self.input_file.read(self.B, "B")
    
    # Thickness    
    self.H = Function(self.V_cg)
    self.input_file.read(self.H, "H")
    
    # Sliding speed
    self.u_b = Function(self.V_cg)
    self.input_file.read(self.u_b, "u_b_0")
    
    # Sliding speed
    self.h = Function(self.V_cg)
    # Potential
    self.phi = Function(self.V_cg)
    # Effective pressure
    self.N = Function(self.V_cg)
    # Flux
    self.q = Function(self.V_cg)
    # Melt rate
    self.m = Function(self.V_cg)
    
    # Potential at 0 pressure
    phi_m = project(pcs['rho_w'] * pcs['g'] * self.B, self.V_cg)
    # Ice overburden pressure
    p_i = project(pcs['rho_i'] * pcs['g'] * self.H, self.V_cg)
    # Potential at overburden pressure
    self.phi0 = project(phi_m + p_i, self.V_cg)
    
    # Vertex coordinates
    self.coords = self.V_cg.tabulate_dof_coordinates().reshape(self.V_cg.dim(), 2)
    self.coords_x = self.coords[:,0]
    self.coords_y = self.coords[:,1]
    
    
    # Get ith time step
  def get_t(self, i):
    if i < self.num_steps:
      attr = self.input_file.attributes("h/vector_" + str(i))
      return attr['timestamp']


  # Get phi at the ith time step
  def get_phi(self, i):
    if i < self.num_steps:
      self.input_file.read(self.phi, "phi/vector_" + str(i))
      return self.phi
      
      
  # Get N at the ith time step
  def get_N(self, i):
    if i < self.num_steps:
      self.input_file.read(self.N, "N/vector_" + str(i))
      return self.N
        
  
  # Get h at the ith time step
  def get_h(self, i):
    if i < self.num_steps:
      self.input_file.read(self.h, "h/vector_" + str(i))
      return self.h
        
  
  # Get m at the i-th time step
  def get_m(self, i):
    if i < self.num_steps:
      self.input_file.read(self.m, "m/vector_" + str(i))
      return self.m
        
        
  # Get the total melt input at the i-th time step
  def get_total_m(self, i):
    if i < self.num_steps:
      self.get_m(i)
      return assemble(self.m * dx)
        
        
  # Compute the volume of water in the sheet at the i-th time
  def get_sheet_volume(self, i):
    if i < self.num_steps:
      self.get_h(i)
      return assemble(self.h * dx)
    
  
  # Write a netcdf file with the results
  def write_netcdf(self, out_file, title):
    root = Dataset(out_file + '.nc', 'w')
    
    ## Dimensions
    
    # Time
    root.createDimension('time', 1)
    # Spatial dimension of model 
    root.createDimension('dim', 2)
    # Number of nodes in mesh
    dim = root.createDimension('index1', self.V_cg.dim())
    
    ## Variables
    
    # Time
    times = root.createVariable('time', 'f8', ('time',))
    times.units = 's'
    times.long_name = 'time'
    #times[0] = 1.0
    
    # Node coordinates
    coords1 = root.createVariable('coords1', 'f8', ('dim', 'index1'))
    coords1.units = 'm'
    coords1.long_name = 'node coordinates'
    coords1[:] = self.coords.T
    
    # Bed
    B = root.createVariable('B', 'f8', ('index1',))
    B.units = 'm'
    B.long_name = 'bed elevation'
    B[:] = self.B.vector()[:]    
    
    # Ice thickness
    H = root.createVariable('H', 'f8', ('index1',))
    H.units = 'm'
    H.long_name = 'ice thickness'
    H[:] = self.H.vector()[:]
    
    # Effective pressure
    N = root.createVariable('N', 'f8', ('time', 'index1',))
    N.units = 'Pa'
    N.long_name = 'effective pressure'
    N[:] = self.N.vector().array().reshape(1, self.V_cg.dim())
    
    # Sheet thickness
    h = root.createVariable('h', 'f8', ('time', 'index1',))
    h.units = 'm'
    h.long_name = 'water sheet thickness'
    h[:] = self.h.vector().array().reshape(1, self.V_cg.dim())
    
    # Water sheet discharge
    q = root.createVariable('q', 'f8', ('time', 'index1',))
    q.units = 'm^2/s'
    q.long_name = 'water sheet discharge'
    q[:] = self.q.vector().array().reshape(1, self.V_cg.dim())
    
    ## Global attributes
    root.title = title
    root.meshtype = 'unstructured'
    root.institution = 'Jacob Z Downs, UM'
    root.source = 'https://github.com/JacobDowns/SheetModel/tree/shmip'
    root.references = 'Schoof et al. 2012, DOI: 10.1017/jfm.2012.165'
    
    root.close() 
    