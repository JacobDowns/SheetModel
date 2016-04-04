from dolfin import *
from dolfin import MPI, mpi_comm_world
from constants import *
from model import *
from phi_solver import *
from h_solver import *

""" Wrapper class for Schoof's constrained sheet model."""

class SheetModel(Model):

  def __init__(self, model_inputs):
    Model.__init__(self, model_inputs)


    ### Initialize model variables
  
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)
    # Cavity height variable
    self.h = Function(self.V_cg)
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
    # Minimum potential constraint for optimization
    self.phi_min = Function(self.V_cg)
    # Maximum potential constraint for optimization
    self.phi_max = Function(self.V_cg)
    # Hydraulic Potential
    self.phi = Function(self.V_cg)
    # Effective pressure
    self.N = Function(self.V_cg)
    # Water pressure
    self.p_w = Function(self.V_cg)
    # Pressure as a fraction of overburden
    self.pfo = Function(self.V_cg)
    
    self.init_model()
    
    
    ### Derive variables we'll need later and setup boundary conditions
    
    # Potential at 0 pressure
    self.phi_m = project(pcs['rho_w'] * pcs['g'] * self.B, self.V_cg)
    # Ice overburden pressure
    self.p_i = project(pcs['rho_i'] * pcs['g'] * self.H, self.V_cg)
    # Potential at overburden pressure
    self.phi_0 = project(self.phi_m + self.p_i, self.V_cg)
    
    
    ### Setup boundary conditions
    
    # If there are boundary conditions specified, use them. Otherwise apply
    # default bc of 0 pressure on the margin
    if 'd_bcs' in self.model_inputs:
      # Dirichlet boundary conditions
      self.d_bcs = self.model_inputs['d_bcs']    
    else :
      # By default, a marker 0f 1 denotes the margin
      self.d_bcs = [DirichletBC(self.V_cg, self.phi_m, self.boundaries, 1)]
      
      
    ### PDE solver parameters
      
    # If the Newton parameters are specified use them. Otherwise use some
    # defaults
    if 'newton_params' in self.model_inputs :
      self.newton_params = self.model_inputs['newton_params']
    else :
      prm = NonlinearVariationalSolver.default_parameters()
      prm['newton_solver']['relaxation_parameter'] = 1.0
      prm['newton_solver']['relative_tolerance'] = 1e-4
      prm['newton_solver']['absolute_tolerance'] = 1e-4
      prm['newton_solver']['error_on_nonconvergence'] = False
      prm['newton_solver']['maximum_iterations'] = 30
      
      self.newton_params = prm
      
      
    ### Setup some stuff necessary for the bfgs optimization   
      
    # If the minimum potential is specified, then use it. Normally one should 
    # probably not do this
    self.phi_min.assign(self.phi_m)
    if 'phi_min' in self.model_inputs:
      self.phi_min.assign(self.model_inputs['phi_min'])
    
    # If the maximum potential is specified, then use it. Normally one should 
    # probably not do this
    self.phi_max.assign(self.phi_0)
    if 'phi_max' in self.model_inputs:
      self.phi_max.assign(self.model_inputs['phi_max'])
      
     # Apply any Dirichlet bcs to the phi_min and phi_max functions
    for bc in self.d_bcs:
      bc.apply(self.phi_min.vector())
      bc.apply(self.phi_max.vector())
    
    # If optimization parameters are specified use them. Otherwise use some 
    # sensible defaults
    if 'opt_params' in self.model_inputs:
      self.opt_params = self.model_inputs['opt_params']
    else:
      self.opt_params = {'tol' : 1e-2, 'scale' : 15}
      
    
    ### Create objects that solve the model equations
    
    # Potential solver
    self.phi_solver = PhiSolver(self)
    # Sheet height solver
    self.h_solver = HSolver(self)
      

    ### Output files
    
    self.h_out = File(self.out_dir + "h.pvd")
    self.phi_out = File(self.out_dir + "phi.pvd")
    self.pfo_out = File(self.out_dir + "pfo.pvd")
    self.m_out = File(self.out_dir + "m.pvd")
    self.u_b_out = File(self.out_dir + "u_b.pvd")
    self.k_out = File(self.out_dir + "k.pvd")
    
    
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
    
    # Write the initial conditions. These can be used as defaults to initialize
    # the model if a simulation crashes and we want to start it again
    self.output_file.write(self.h, "h_0")
    self.output_file.write(self.m, "m_0")
    self.output_file.write(self.u_b, "u_b_0")
    self.output_file.write(self.mesh, "mesh")
    self.output_file.write(self.B, "B")
    self.output_file.write(self.H, "H")
    self.output_file.write(self.k, "k_0")
    self.output_file.write(self.boundaries, "boundaries")
    
    # Initial time is 0
    self.t = 0.0        
    
    
  # Initialize model using the state at the end of a previous simulation
  def continue_simulation(self):
     self.load_inputs()
     
     # Load the most recent cavity height and channel height values
     num_steps = self.input_file.attributes("h")['count']
     h_last = "h/vector_" + str(num_steps - 1)
     self.input_file.read(self.h, h_last)
      
     # Get the start time for the simulation        
     attr = self.input_file.attributes(h_last)
     self.t = attr['timestamp']
    
    
  # Steps phi and h forward by dt
  def step(self, dt):
    self.phi_solver.step()
    self.h_solver.step(dt)
    # Update model time
    self.t += dt
    
  
  # Steps phi forward using the optimization procedure then steps h 
  # forward
  def step_opt(self, dt):
    self.phi_solver.solve_opt()
    self.h_solver.step(dt)
    # Update model time
    self.t += dt
    
    
  # Load all model inputs from the input file if they are not specified
  # in the
  def load_inputs(self):
    try :
      # Bed elevation 
      self.input_file.read(self.B, "B")
      # Ice thickness
      self.input_file.read(self.H, "H")    
      # Boundary facet function
      self.input_file.read(self.boundaries, "boundaries")
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
      self.k.assign(project(self.k, self.V_cg))
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
    self.update_pfo()
    
  
  # Update the pressure as a fraction of overburden to reflect the current 
  # value of phi
  def update_pfo(self):
    # Compute overburden pressure
    self.pfo.vector().set_local(self.p_w.vector().array() / self.p_i.vector().array())
    self.pfo.vector().apply("insert")
    
  
  # Updates all fields derived from phi
  def update_phi(self):
    self.update_N()
    self.update_pw()
    
  
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
      if 'pfo' in to_write:
        self.pfo_out << self.pfo
      if 'u' in to_write:
        self.u_out << self.phi_solver.u
      if 'm' in to_write:
        self.m_out << self.m
      if 'u_b' in to_write:
        self.u_b_out << self.u_b
      if 'q' in to_write:
        self.q_func.assign(project(self.phi_solver.q, self.V_cg))
        self.q_out << self.q_func
      if 'k' in to_write:
        self.k_out << self.k


  # Write checkpoint files to an hdf5 file
  def checkpoint(self, to_write = [], label = None):
    to_write = set(to_write)

    # Always checkpoint h
    self.output_file.write(self.h, "h", self.t)

    if 'phi' in to_write:
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
    
  
  # Updates the melt rate function
  def set_m(self, new_m):
    self.m.assign(new_m)
    self.m_func.assign(new_m)
    
  
  # Update sliding speed
  def set_u_b(self, new_u_b):
    self.u_b.assign(new_u_b)
    self.u_b_func.assign(new_u_b)
    
    
  # Sets the hydraulic conductivity
  def set_k(self, new_k):
    self.k.assign(project(new_k, self.V_cg))
    self.update_k()
    