from dolfin import *
from dolfin import MPI, mpi_comm_world
from cr_tools import *
from model import *
from phi_solver import *
from hs_solver import *
import sys as sys

""" Wrapper class for Glads."""

class ChannelModel(Model):

  def __init__(self, model_inputs):
    Model.__init__(self, model_inputs)

    ### Initialize model variables

    # CG function space
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)
    # CR function space
    self.V_cr = FunctionSpace(self.mesh, "CR", 1)
    # Cavity height
    self.h = Function(self.V_cg)
    # Chanel height
    self.S = Function(self.V_cr)
    # Bed geometry
    self.B = Function(self.V_cg)
    # Ice thickness
    self.H = Function(self.V_cg)
    # Melt input
    self.m = Function(self.V_cg)
    # Basal sliding speed
    self.u_b = Function(self.V_cg)
    # Hydraulic conductivity
    self.k = Function(self.V_cg)
    # Bump height 
    self.h_r = Function(self.V_cg)
    # Boundary facet function
    self.boundaries = FacetFunction("size_t", self.mesh)
    # Hydraulic Potential
    self.phi = Function(self.V_cg)
    # Potential at previous time step
    self.phi_prev = Function(self.V_cg)
    # Effective pressure
    self.N = Function(self.V_cg)
    # Water pressure
    self.p_w = Function(self.V_cg)
    # Pressure as a fraction of overburden
    self.pfo = Function(self.V_cg)
    # A cr function mask that is 1 on interior edges and 0 on exterior edges. 
    # Used to prevent opening on exterior edges.
    self.mask = Function(self.V_cr)
    # Derivative of potential over channel edges
    self.dphi_ds_cr = Function(self.V_cr)
    # Effective pressure on edges
    self.N_cr = Function(self.V_cr)
    # Sheet height on edges
    self.h_cr = Function(self.V_cr)
    # Hydraulic conductivity on edges
    self.k_cr = Function(self.V_cr)
    # Stores the value of S**alpha. A workaround for a bug in Fenics that
    # causes problems when exponentiating a CR function
    self.S_alpha = Function(self.V_cr)
    # Edge lengths 
    self.edge_lens = Function(self.V_cr)


    # Load inputs and initialize input and output files    
    self.init_model()
    # CR tools for dealing with CR functions
    self.cr_tools = CRTools(self.mesh, self.V_cg, self.V_cr, self.edge_lens)
    
    ### Derive some additional fields
    
    # Potential at 0 pressure
    self.phi_m = project(pcs['rho_w'] * pcs['g'] * self.B, self.V_cg)
    # Ice overburden pressure
    self.p_i = project(pcs['rho_i'] * pcs['g'] * self.H, self.V_cg)
    # Potential at overburden pressure
    self.phi_0 = project(self.phi_m + self.p_i, self.V_cg)
    
    # Populate all fields derived from the primary variables
    self.update_phi()
    self.update_h()
    self.update_S()
    self.update_k()


    ### Setup boundary conditions
    
    # If there are boundary conditions specified, use them. Otherwise apply
    # default bc of 0 pressure on the margin
    if 'd_bcs' in self.model_inputs:
      # Dirichlet boundary conditions
      self.d_bcs = self.model_inputs['d_bcs']    
    else :
      # By default, a marker 0f 1 denotes the margin
      self.d_bcs = [DirichletBC(self.V_cg, self.phi_m, self.boundaries, 1)]
      
      
    ### Newton solver parameters
      
    # If the Newton parameters are specified use them. Otherwise use some
    # defaults
    if 'newton_params' in self.model_inputs :
      self.newton_params = self.model_inputs['newton_params']
    else :
      prm = NonlinearVariationalSolver.default_parameters()
      prm['newton_solver']['relaxation_parameter'] = 1.0
      prm['newton_solver']['relative_tolerance'] = 1e-3
      prm['newton_solver']['absolute_tolerance'] = 1e-3
      prm['newton_solver']['error_on_nonconvergence'] = False
      prm['newton_solver']['maximum_iterations'] = 25
      
      self.newton_params = prm
      
      
    ### Create objects that solve the model equations
    
    # Potential solver
    self.phi_solver = PhiSolver(self)
    # Sheet height solver
    self.hs_solver = HSSolver(self)
      
      
    ### Create output files
    
    # Output directory
    self.h_out = File(self.out_dir + "h.pvd")
    self.S_out = File(self.out_dir + "S.pvd")
    self.phi_out = File(self.out_dir + "phi.pvd")
    self.N_out = File(self.out_dir + "N.pvd")
    self.pfo_out = File(self.out_dir + "pfo.pvd")
    self.m_out = File(self.out_dir + "m.pvd")
    self.u_b_out = File(self.out_dir + "u_b.pvd")
    self.k_out = File(self.out_dir + "k.pvd")
    self.h_cr_out = File(self.out_dir + "h_cr.pvd")    
    self.N_cr_out = File(self.out_dir + "N_cr.pvd")
    
    # Facet functions for plotting CR functions in Paraview    
    self.ff_out_S = FacetFunctionDouble(self.mesh)
    self.ff_out_N_cr = FacetFunctionDouble(self.mesh)
    self.ff_out_h_cr = FacetFunctionDouble(self.mesh)  
    
    
  # Steps phi and h forward by dt
  def step(self, dt):
    self.phi_solver.step(dt)
    self.hs_solver.step(dt)

    # Update model time
    self.t += dt
    
  
  # Look at the input file to check if we're starting or continuing a simulation
  def to_continue(self):
    try :
      # Check if this is a simulation output file. If it is it should have an 
      # attribute called h
      self.input_file.attributes("h")
      return True
    except :
      return False
  
    
  # Initialize the model state for a new simulation
  def start_simulation(self):
    self.load_inputs()
  
    # Read in the initial cavity height
    self.input_file.read(self.h, "h_0")
    # Load the intitial channel heights
    self.input_file.read(self.S, "S_0")
    # Load the initial hydraulic potential
    self.input_file.read(self.phi, "phi_0")
    
    # Write the initial conditions. These can be used as defaults to initialize
    # the model if a simulation crashes and we want to start it again
    self.output_file.write(self.h, "h_0")
    self.output_file.write(self.m, "m_0")
    self.output_file.write(self.u_b, "u_b_0")
    self.output_file.write(self.phi, "phi_0")
    self.output_file.write(self.mask, "mask")
    self.output_file.write(self.edge_lens, "edge_lens")
    self.output_file.write(self.mesh, "mesh")
    self.output_file.write(self.B, "B")
    self.output_file.write(self.H, "H")
    self.output_file.write(self.boundaries, "boundaries")
    self.output_file.write(self.S, "S_0")    
    
  
  # Initialize model using the state at the end of a previous simulation
  def continue_simulation(self):
     self.load_inputs()
     
     # Load the most recent cavity height and channel height values
     num_steps = self.input_file.attributes("h")['count']
     h_last = "h/vector_" + str(num_steps - 1)
     S_last = "S/vector_" + str(num_steps - 1)
     phi_last = "phi/vector_" + str(num_steps - 1)
      
     self.input_file.read(self.h, h_last)
     self.input_file.read(self.S, S_last)
     self.input_file.read(self.phi, phi_last)
      
     # Get the start time for the simulation        
     attr = self.input_file.attributes(h_last)
     self.t = attr['timestamp']
     
     
  # Load all model inputs from the input file if they are not specified
  # in the
  def load_inputs(self):
    try :
      # Bed elevation 
      self.input_file.read(self.B, "B")
      # Ice thickness
      self.input_file.read(self.H, "H")    
      # Edge mask
      self.input_file.read(self.mask, "mask")
      # Edge lengths
      self.input_file.read(self.edge_lens, "edge_lens")
      # Boundary facet function
      self.input_file.read(self.boundaries, "boundaries")
      # Hydraulic potential
      self.assign_func(self.phi, "phi_0")
      self.phi_prev.assign(self.phi)
      # Melt input
      self.assign_func(self.m, "m_0")
      # Sliding speed
      self.assign_func(self.u_b, "u_b_0")
  
    except Exception as e:
      # If we can't find one of these model inputs we're done 
      print >> sys.stderr, e
      print >> sys.stderr, "Could not load model inputs."
      sys.exit(1)
    
    try :
      # If we get a hydraulic conductivity expression or function use it
      self.assign_func(self.k, 'k')
    except :
      # Otherwise we'll just use use the default constant conductivity 
      self.k.assign(interpolate(Constant(self.pcs['k']), self.V_cg))
      
    try :
      # If we get a bump height function use it
      self.assign_func(self.h_r, 'h_r')
    except :
      # Otherwise we'll just use use the default constant bump height
      self.h_r.assign(interpolate(Constant(self.pcs['h_r']), self.V_cg))
      
      
  # Update the effective pressure to reflect current value of phi
  def update_N(self):
    self.N.vector().set_local(self.phi_0.vector().array() - self.phi.vector().array())
    self.N.vector().apply("insert")
    
  
  # Update the water pressure to reflect current value of phi
  def update_pw(self):
    self.p_w.vector().set_local(self.phi.vector().array() - self.phi_m.vector().array())
    self.p_w.vector().apply("insert")
    
  
  # Update the pressure as a fraction of overburden to reflect the current 
  # value of phi
  def update_pfo(self):
    # Update water pressure
    self.update_pw()
  
    # Compute overburden pressure
    self.pfo.vector().set_local(self.p_w.vector().array() / self.p_i.vector().array())
    self.pfo.vector().apply("insert")
    
  
  # Updates all fields derived from phi
  def update_phi(self):
    self.update_N()
    self.update_pfo()
    self.update_dphi_ds_cr()
    self.update_N_cr()


  # Updates all fields derived from h    
  def update_h(self):
    self.update_h_cr()
    
    
  # Updates all fields derived from S
  def update_S(self):
    self.update_S_alpha()
    
    
  # Update the edge derivatives of the potential to reflect current value of phi
  def update_dphi_ds_cr(self):
    self.cr_tools.ds(self.phi, self.dphi_ds_cr)
    
  
  # Update effective pressure on edge midpoints to reflect current value of phi
  def update_N_cr(self):
    self.update_N()
    self.cr_tools.midpoint(self.N, self.N_cr)
    
  
  # Update the edge midpoint values h_cr to reflect the current value of h
  def update_h_cr(self):
    self.cr_tools.midpoint(self.h, self.h_cr)
    
  
  # Update S**alpha to reflect current value of S
  def update_S_alpha(self):
    alpha = self.pcs['alpha']
    self.S_alpha.vector().set_local(self.S.vector().array()**alpha)
    
  # Update fields derived from k
  def update_k(self):
    self.cr_tools.midpoint(self.k, self.k_cr)
    
  
  # Write fields to pvd files for visualization
  def write_pvds(self, to_write = []):
    to_write = set(to_write)
    if len(to_write) == 0:
      self.h_out << self.h
      self.pfo_out << self.pfo
    else:
      if 'h' in to_write:
        self.h_out << self.h
      if 'phi' in to_write:
        self.h_out << self.h
      if 'N' in to_write:
        self.N_out << self.N
      if 'pfo' in to_write:
        self.pfo_out << self.pfo
      if 'u' in to_write:
        self.u_out << self.phi_solver.u
      if 'm' in to_write:
        self.m_out << self.m_func
      if 'u_b' in to_write:
        self.u_b_out << self.u_b_func
      if 'q' in to_write:
        self.q_func.assign(project(self.phi_solver.q, self.V_cg))
        self.q_out << self.q_func
      if 'k' in to_write:
        self.k_out << self.k_func
      if 'S' in to_write:
        self.cr_tools.copy_cr_to_facet(self.S, self.ff_out_S)
        self.S_out << self.ff_out_S
      if 'h_cr' in to_write:
        self.cr_tools.copy_cr_to_facet(self.h_cr, self.ff_out_h_cr)
        self.h_cr_out << self.ff_out_h_cr
      if 'N_cr' in to_write:
        self.cr_tools.copy_cr_to_facet(self.N_cr, self.ff_out_N_cr)
        self.N_cr_out << self.ff_out_N_cr
      

  # Write checkpoint files to an hdf5 file
  def checkpoint(self, to_write = [], label = None):
    to_write = set(to_write)

    # Always checkpoint h and S
    self.output_file.write(self.h, "h", self.t)
    self.output_file.write(self.S, "S", self.t)
    self.output_file.write(self.phi, "phi", self.t)

    if 'u_b' in to_write:
      self.output_file.write(self.u_b, "u_b", self.t)
    if 'k' in to_write:
      self.output_file.write(self.k, "k", self.t)
    if 'm' in to_write:
      self.output_file.write(self.m, "m", self.t)
    if 'pfo' in to_write:
      self.output_file.write(self.pfo, "pfo", self.t)
      
    self.output_file.flush()
    
  
  # Sets the melt rate function
  def set_m(self, new_m):
    self.m.assign(new_m)
    
  
  # Sets the sliding speed
  def set_u_b(self, new_u_b):
    self.u_b.assign(new_u_b)
    
    
  # Sets the hydraulic conductivity
  def set_k(self, new_k):
    self.k.assign(new_k)
    self.update_k()
    
    