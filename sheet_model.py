from dolfin import *
from dolfin_adjoint import *
from constants import *
from phi_solver import *
from h_solver import *

""" Wrapper class for Schoof's constrained sheet model."""

class SheetModel():

  def __init__(self, model_inputs):

    ### Initialize model inputs

    # Dictionary of model inputs
    self.model_inputs = model_inputs
    
    # If the model is initialized in start mode, then create a new output file
    if not 'mode' in model_inputs or model_inputs['mode'] == 'start' :

      # Get the output file name if there is one
      out_name = 'out'
      if 'output_file' in model_inputs:
        out_name = model_inputs['output_file']
        
      # Load the input file 
      self.input_file = HDF5File(mesh.mpi_comm(), self.model_inputs['input_file'], 'r')
      # Create an output file 
      self.ouput_file = HDF5File(mesh.mpi_comm(), out_name + ".hdf5", 'w')
    else :
      # If not initialized in start mode, then we assume that the user wants to 
      # continue an existing simulation
    
      # Load the input file in append mode
      self.input_file = HDF5File(mesh.mpi_comm(), self.model_inputs['input_file'], 'a')
      # The output file will be the same
      self.ouput_file = self.input_file
    
    # Load the mesh
    self.mesh = self.input_file.read(mesh, "mesh", False)    
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)

    # If there is a dictionary of physical constants specified in model_inputs, 
    # use it. Otherwise use the defaults. 
    if 'constants' in self.model_inputs :
      self.pcs = self.model_inputs['constants']
    else :
      self.pcs = pcs
    
    # If an input directory is specified, load model inputs from there. 
    # Otherwise use the specified model inputs dictionary.
    if in_dir:
      self.load_inputs(in_dir)
      
      
    # Ice thickness    
    self.H = Function(self.V_cg)
    self.io_file.read(self.H, "H")
    
    # Bed elevation     
    self.B = Function(self.V_cg)
    self.io_file.read(self.B, "B")
    
    # Basal sliding speed
    self.u_b = Function(self.V_cg)
    self.u_b_func = Function(self.V_cg)
    self.assign_var('u_b', self.u_b, self.u_b_func)
    
    # Melt rate
    self.m = Function(self.V_cg)
    self.m_func = Function(self.V_cg)
    self.assign_var('m', self.m, self.m_func)
    
    # Cavity gap height
    self.h = Function(self.V_cg)
    self.io_file.
    # Potential at 0 pressure
    self.phi_m = project(pcs['rho_w'] * pcs['g'] * self.B, self.V_cg)
    # Ice overburden pressure
    self.p_i = project(pcs['rho_i'] * pcs['g'] * self.H, self.V_cg)
    # Potential at overburden pressure
    self.phi_0 = project(self.phi_m + self.p_i, self.V_cg)


    # Facet function marking boundaries
    self.boundaries = self.model_inputs['boundaries']
    
    # If there are boundary conditions specified, use them. Otherwise apply
    # default bc of 0 pressure on the margin
    if 'd_bcs' in self.model_inputs:
      # Dirichlet boundary conditions
      self.d_bcs = self.model_inputs['d_bcs']    
    else :
      # By default, a marker 0f 1 denotes the margin
      self.d_bcs = [DirichletBC(self.V_cg, self.phi_m, self.boundaries, 1)]

    # Hydraulic conductivity
    self.k = Function(self.V_cg)
    self.k_func = Function(self.V_cg)
        
    # If we get a hydraulic conductivity expression or function use it
    if 'k' in self.model_inputs:
      self.assign_var('k', self.k, self.k_func)
    else :
      # Otherwise use the default constant conductivity 
      self.k = Function(self.V_cg)
      self.k.interpolate(Constant(self.pcs['k']))
      self.k_func.assign(self.k)
      
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

    # If the initial time is specified then use it, otherwise use the default
    # initial time of 0

    # Current time
    self.t = 0.0
    if 't0' in self.model_inputs:
      self.t = self.model_inputs['t0']


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
    self.out_dir = self.model_inputs['out_dir']
    self.h_out = File(self.out_dir + "h.pvd")
    self.phi_out = File(self.out_dir + "phi.pvd")
    self.pfo_out = File(self.out_dir + "pfo.pvd")
    self.u_out = File(self.out_dir + "u.pvd")
    self.m_out = File(self.out_dir + "m.pvd")
    self.u_b_out = File(self.out_dir + "u_b.pvd")
    self.q_out = File(self.out_dir + "q.pvd")
    self.k_out = File(self.out_dir + "k.pvd")
    # An index for labeling checkpoint files
    self.check_index = 0
    # Directory to write checkpoint files
    self.check_dir = self.out_dir + "checkpoint/"


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
    
    
  # Load all model inputs from a directory except for the mesh and initial 
  # conditions on h, h_w, and phi
  def load_inputs(self, in_dir):
      
    # Melt
    if not 'm' in self.model_inputs:
      m = Function(self.V_cg)
      File(in_dir + "m.xml") >> m
      self.model_inputs['m'] = m
      
    # Sliding speed
    if not 'u_b' in self.model_inputs:
      u_b = Function(self.V_cg)
      File(in_dir + "u_b.xml") >> u_b 
      self.model_inputs['u_b'] = u_b
      
    # Boundary facet function
    if not 'boundaries' in self.model_inputs:
      boundaries = FacetFunction('size_t', self.mesh)
      File(in_dir + "boundaries.xml") >> boundaries
      self.model_inputs['boundaries'] = boundaries
    
  
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
    
  
  # Write fields to pvd files
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


  # Write out xml checkpoint flies for h and phi. to_write is a list of the fields
  # to output. label is simply a way to keep track of which file is which (could
  # be an index or a time or something)
  def write_xmls(self, to_write = [], label = None):
    to_write = set(to_write)
    
    if not label:
      label = self.check_index
    
    if 'h' in to_write:
      File(self.check_dir + "h_" + str(label) + ".xml") << self.h
    if 'phi' in to_write:
      File(self.check_dir + "phi_" + str(label) + ".xml") << self.phi
    if 'u' in to_write:
      File(self.check_dir + "u_" + str(label) + ".xml") << self.phi_solver.u
    if 'u_b' in to_write:
      File(self.check_dir + "u_b_" + str(label) + ".xml") << self.u_b_func
    if 'k' in to_write:
      File(self.check_dir + "k_" + str(label) + ".xml") << self.k_func
    if 'pfo' in to_write:
      File(self.check_dir + "pfo_" + str(label) + ".xml") << self.pfo
        
    # Increment the checkpoint index
    self.check_index += 1
    
  
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
    
  
  # Some variables need to be functions for writing to a file or other reasons.
  # This function takes in a variable that could be either a function or expression
  # and assigns it to a function form of the variable (var_func)
  def assign_var_func(self, var, var_func):
    if isinstance(var, dolfin.Expression):
        var_func.assign(project(var, self.V_cg))
    else :
        var_func.assign(var)
   
   
  # Assigns a variable that could be an expression or a function. The function
  # first checks the model_inputs dictionary, then the hdf5 io file for the variable 
  def assign_var(self, name, var, var_func):
    # Look for the variable in the model_inputs directory
    if 'name' in self.model_inputs:
      input_var = self.model_inputs[name]
    else:
      # Otherwise get it from the io_file
      input_var = self.io_file[name]
    
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