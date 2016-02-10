from dolfin import *
from dolfin import MPI, mpi_comm_world
from dolfin_adjoint import *
from constants import *
from phi_solver import *
from h_solver import *

""" Wrapper class for Schoof's constrained sheet model."""

class SheetModel():

  def __init__(self, model_inputs):

    ### Load some model inputs

    # Dictionary of model inputs
    self.model_inputs = model_inputs

    # Check if an input file is specified . This input_file should minimally
    # contain bed elevation, ice thickness, and a facet function for marking
    # boundaries. It can contain all model inputs if desired.
    if 'input_file' in self.model_inputs:
      self.input_file = HDF5File(mpi_comm_world(), self.model_inputs['input_file'], 'r') 
    else :
      # Otherwise yell at the user (that would be myself)
      raise Exception("Please specify a valid input file, if you would be so kind.")
      sys.exit(1)
      
    # Output directory for writing xdmf files for visualization and the checkpoint
    # file, if we're not continuing a simulation
    if 'out_dir' in self.model_inputs:
      self.out_dir = model_inputs['out_dir'] + '/'
    else :
      # Otherwise yell at the user (that would be myself)
      raise Exception("Please specify an output directory, if you would be so kind.")
      sys.exit(1)

    # Load the mesh
    self.mesh = Mesh()
    self.input_file.read(self.mesh, "mesh", False)    
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)

    # Cavity height variable
    self.h = Function(self.V_cg)

    # Bed geometry
    self.B = Function(self.V_cg)

    # Ice thickness
    self.H = Function(self.V_cg)
    
    # Melt input
    self.m = Function(self.V_cg)
    # Function form of melt
    self.m_func = Function(self.V_cg)
    
    # Basal sliding speed
    self.u_b = Function(self.V_cg)
    # Function form of sliding speed
    self.u_b_func = Function(self.V_cg)
    
    # Hydraulic conductivity
    self.k = Function(self.V_cg)
    # Function form of conductivity
    self.k_func = Function(self.V_cg)
    
    # Boundary facet function
    self.boundaries = FacetFunction("size_t", self.mesh)
    
    # If there is a dictionary of physical constants specified in model_inputs, 
    # use it. Otherwise use the defaults. 
    if 'constants' in self.model_inputs :
      self.pcs = self.model_inputs['constants']
    else :
      self.pcs = pcs
        
    # Load a bunch of model inputs including the ice geometry, melt, sliding, and 
    # conductivity
    self.load_inputs()
    
    
    ### Derive variables we'll need later and setup boundary conditions
    
    # Potential at 0 pressure
    self.phi_m = project(pcs['rho_w'] * pcs['g'] * self.B, self.V_cg)
    # Ice overburden pressure
    self.p_i = project(pcs['rho_i'] * pcs['g'] * self.H, self.V_cg)
    # Potential at overburden pressure
    self.phi_0 = project(self.phi_m + self.p_i, self.V_cg)
    
    # If there are boundary conditions specified, use them. Otherwise apply
    # default bc of 0 pressure on the margin
    if 'd_bcs' in self.model_inputs:
      # Dirichlet boundary conditions
      self.d_bcs = self.model_inputs['d_bcs']    
    else :
      # By default, a marker 0f 1 denotes the margin
      self.d_bcs = [DirichletBC(self.V_cg, self.phi_m, self.boundaries, 1)]
      
      
    ### Numerical solver parameters
      
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
      prm['newton_solver']['maximum_iterations'] = 25
      
      self.newton_params = prm
      
    # If optimization parameters are specified use them. Otherwise use some 
    # sensible defaults
    if 'opt_params' in self.model_inputs:
      self.opt_params = self.model_inputs['opt_params']
    else:
      self.opt_params = {'tol' : 0.5e-7, 'scale' : 15.0}
      
      
    # Figure out if we should start a new simulation or continue running an 
    #  existing one based on what the input file looks like 
    continue_simulation = False
        
    try :
      # Check if this is a simulation output file. If it is it should have an 
      # attribute called h. If this isn't the case, this will throw an exception
      self.input_file.attributes("h")
      continue_simulation = True
    except :
      pass
    
    if continue_simulation : 
      # Continue an existing simulation
      
      # Close the input file and re-open in append mode
      self.input_file.close()
      self.input_file = HDF5File(mpi_comm_world(), self.model_inputs['input_file'], 'a') 
      
      # Load the most recent cavity height value
      num_steps = self.input_file.attributes("h")['count']
      h_last = "h/vector_" + str(num_steps - 1)
      self.input_file.read(self.h, h_last)
      
      # Get the end time for the simulation        
      attr = self.input_file.attributes(h_last)
      self.t = attr['timestamp']
      
      # We'll just checkpoint by appending the existing output file rather than 
      # making a new one 
      self.output_file = self.input_file
    else :
      # Start a new simulation 
    
      # Read in the initial condition
      self.input_file.read(self.h, "h_0")
      
      # Check if there's a given name for the checkpoint file
      checkpoint_file = 'out'
      if 'checkpoint_file' in model_inputs:
        checkpoint_file = model_inputs['checkpoint_file']
      
      # Create a new ouput file for checkpointing.
      self.output_file = HDF5File(mpi_comm_world(), self.out_dir + checkpoint_file + ".hdf5", 'w')
      
      # Write the initial conditions. These can be used as defaults to initialize
      # the model if a simulation crashes and we want to start it again
      self.output_file.write(self.h, "h_0")
      self.output_file.write(self.m, "m_0")
      self.output_file.write(self.u_b, "u_b_0")
      self.output_file.write(self.mesh, "mesh")
      # And the ice sheet geometry
      self.output_file.write(self.B, "B")
      self.output_file.write(self.H, "H")
      self.output_file.write(self.boundaries, "boundaries")
      
      # Initial time is 0
      self.t = 0.0        
        

    ### Set up a few more things we'll need

    # Potential
    self.phi = Function(self.V_cg)
    # Effective pressure
    self.N = Function(self.V_cg)
    # Water pressure
    self.p_w = Function(self.V_cg)
    # Pressure as a fraction of overburden
    self.pfo = Function(self.V_cg)
    
    
    ### Output files
    
    # Output directory
    self.h_out = File(self.out_dir + "h.pvd")
    self.phi_out = File(self.out_dir + "phi.pvd")
    self.pfo_out = File(self.out_dir + "pfo.pvd")
    self.m_out = File(self.out_dir + "m.xdmf")
    self.u_b_out = File(self.out_dir + "u_b.pvd")
    self.k_out = File(self.out_dir + "k.pvd")
    

    ### Create the solver objects
    
    # Potential solver
    self.phi_solver = PhiSolver(self)
    # Sheet height solver
    self.h_solver = HSolver(self)
    
    
  # Steps phi and h forward by dt
  def step(self, dt):
    self.update_time_dependent_vars()
    self.phi_solver.step()
    self.h_solver.step(dt)
    # Update model time
    self.t += dt
    
  
  # Steps phi forward using the optimization procedure then steps h 
  # forward
  def step_opt(self, dt):
    self.update_time_dependent_vars()
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
      self.assign_var("m_0", self.m, self.m_func)
      
      # Sliding speed
      self.assign_var("u_b_0", self.u_b, self.u_b_func)     
    except Exception as e:
      # If we can't find one of these model inputs we're done 
      print >> sys.stderr, e
      print >> sys.stderr, "Could not load model inputs."
      sys.exit(1)
    
    try :
      # If we get a hydraulic conductivity expression or function use it
      self.assign_var('k', self.k, self.k_func)
    except :
      # Otherwise we'll just use use the default constant conductivity 
      self.k.interpolate(Constant(self.pcs['k']))
      self.k_func.assign(self.k)
      
      
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
        self.m_out << self.m_func
      if 'u_b' in to_write:
        self.u_b_out << self.u_b_func
      if 'q' in to_write:
        self.q_func.assign(project(self.phi_solver.q, self.V_cg))
        self.q_out << self.q_func
      if 'k' in to_write:
        self.k_out << self.k_func


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
      
    self.output_file.flush()
    
  
  # Updates the melt rate function
  def set_m(self, new_m):
    self.m.assign(new_m)
    self.m_func.assign(new_m)
    
  
  # Update sliding speed
  def set_u_b(self, new_u_b):
    self.u_b.assign(new_u_b)
    self.u_b_func.assign(new_u_b)
    
    
  # Updates the hydraulic conductivity
  def set_k(self, new_k):
    self.k.assign(new_k)
    self.k_func.assign(new_k)
   
   
  # Assigns a value to a model variable. This function first looks in the 
  # model_inputs dictionary then the input_file. 
  def assign_var(self, name, var, var_func):
    # Look for the variable in the model_inputs dictionary. It could be either
    # a function or expression.
    if 'name' in self.model_inputs:
      input_var = self.model_inputs[name]
      
      if isinstance(input_var, dolfin.Expression):
        # If the input variable is an expression, we need to project it to 
        # get a function
        var = input_var
        var_func.assign(project(input_var, self.V_cg))
      else :
        # If the input variable is a function, then we'll copy it rather than
        # use the original function to prevent unwanted metaphysical linkage
        var.assign(input_var)
        var_func.assign(input_var)
          
    else :
      # If we didn't find anything in the model_inputs dictionary, check the input_file.
      if self.read_var(name, var) :
        var_func.assign(var)
      else :
        # If the input isn't there either raise an exception
        raise Exception("Could not load model input: " + str(name))
    
  
  # Reads in a variable from 
  def read_var(self, name, var):
    try :
      self.input_file.read(var, name)
      return True
    except :
      return False
    
  # Updates the potentially time dependent expressions (m, u_b, k)
  def update_time_dependent_vars(self):
    # Check if the melt is an expression and if it has an attribute t. If so 
    # update t to the model time
    if isinstance(self.m, dolfin.Expression):
      if hasattr(self.m, 't'):
        # Update the time
        self.m.t = self.t  
        # Update function form of m
        self.m_func.assign(project(self.m, self.V_cg))
        
    # Same deal for sliding speed
    if isinstance(self.u_b, dolfin.Expression):
      if hasattr(self.u_b, 't'):
        # Update the time
        self.u_b.t = self.t  
        # Update function form of u_b
        self.u_b_func.assign(project(self.u_b, self.V_cg))
        
    # Same deal for k
    if isinstance(self.k, dolfin.Expression):
      if hasattr(self.k, 't'):
        self.k.t = self.t  
        self.k_func.assign(project(self.k, self.V_cg))
