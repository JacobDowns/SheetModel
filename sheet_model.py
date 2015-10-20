from dolfin import *
from dolfin_adjoint import *
from constants import *
from phi_solver import *
from h_solver import *

""" Wrapper class Schoof's constrained sheet model."""

class SheetModel():

  def __init__(self, model_inputs, in_dir = None):

    ### Initialize model inputs

    self.mesh = model_inputs['mesh']
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)
    self.model_inputs = model_inputs
    
    # If an input directory is specified, load model inputs from there. 
    # Otherwise use the specified model inputs dictionary.
    if in_dir:
      self.load_inputs(in_dir)
      
    # Ice thickness    
    self.H = self.model_inputs['H']
    # Bed elevation       
    self.B = self.model_inputs['B']
    # Basal sliding speed
    self.u_b = self.model_inputs['u_b']
    # Function form of sliding speed
    self.u_b_func = Function(self.V_cg)
    if isinstance(self.u_b, dolfin.Expression):
        self.u_b_func.assign(project(self.u_b, self.V_cg))
    else :
        self.u_b_func.assign(self.u_b)
    # Melt rate
    self.m = self.model_inputs['m']
    # Function form of melt rate
    self.m_func = Function(self.V_cg)
    if isinstance(self.m, dolfin.Expression):
        self.m_func.assign(project(self.m, self.V_cg))
    else :
        self.m_func.assign(self.m)
    # Cavity gap height
    self.h = self.model_inputs['h_init']
    # Potential at 0 pressure
    self.phi_m = self.model_inputs['phi_m']
    # Ice overburden pressure
    self.p_i = self.model_inputs['p_i']
    # Potential at overburden pressure
    self.phi_0 = self.model_inputs['phi_0']
    # Dirichlet boundary conditions
    self.d_bcs = self.model_inputs['d_bcs']
    # Facet function marking boundaries
    self.boundaries = self.model_inputs['boundaries']
    # Output directory
    self.out_dir = self.model_inputs['out_dir']
    
    # If there is a dictionary of physical constants specified, use it. 
    # Otherwise use the defaults. 
    if 'constants' in self.model_inputs :
      self.pcs = self.model_inputs['constants']
    else :
      self.pcs = pcs
      
    # If the Newton parameters are specified use them. Otherwise use some
    # defaults
    if 'newton_params' in self.model_inputs :
      self.newton_params = self.model_inputs['newton_params']
    else :
      prm = NonlinearVariationalSolver.default_parameters()
      prm['newton_solver']['relaxation_parameter'] = 1.
      prm['newton_solver']['relative_tolerance'] = 3e-6
      prm['newton_solver']['absolute_tolerance'] = 1e-3
      prm['newton_solver']['error_on_nonconvergence'] = False
      prm['newton_solver']['maximum_iterations'] = 25
      
      self.newton_params = prm
      
    # If optimization parameters are specified use them. Otherwise use some 
    # sensible defaults
    if 'opt_params' in self.model_inputs:
      self.opt_params = self.model_inputs['opt_params']
    else:
      self.opt_params = {'tol' : 2e-8, 'scale' : 150.0}


    ### Set up a few more things we'll need

    # Potential
    self.phi = Function(self.V_cg)
    # Effective pressure
    self.N = Function(self.V_cg)
    # Water pressure
    self.p_w = Function(self.V_cg)
    # Pressure as a fraction of overburden
    self.pfo = Function(self.V_cg)
    # Current time
    self.t = 0.0
    
    
    ### Output files
    
    self.h_out = File(self.out_dir + "h.pvd")
    self.phi_out = File(self.out_dir + "phi.pvd")
    self.pfo_out = File(self.out_dir + "pfo.pvd")
    self.u_out = File(self.out_dir + "u.pvd")
    self.m_out = File(self.out_dir + "m.pvd")
    self.u_b_out = File(self.out_dir + "u_b.pvd")
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
    # Bed
    if not 'B' in self.model_inputs:
      B = Function(self.V_cg)
      File(in_dir + "B.xml") >> B
      self.model_inputs['B'] = B
      
    # Ice thickness
    if not 'H' in self.model_inputs:
      H = Function(self.V_cg)
      File(in_dir + "H.xml") >> H
      self.model_inputs['H'] = H
      
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
      
    # Potential at 0 pressure
    if not 'phi_m' in self.model_inputs:
      phi_m = Function(self.V_cg)
      File(in_dir + "phi_m.xml") >> phi_m
      self.model_inputs['phi_m'] = phi_m
      
    # Potential at overburden pressure
    if not 'phi_0' in self.model_inputs:
      phi_0 = Function(self.V_cg)
      File(in_dir + "phi_0.xml") >> phi_0
      self.model_inputs['phi_0'] = phi_0
      
    # Ice overburden pressure
    if not 'p_i' in self.model_inputs:
      p_i = Function(self.V_cg)
      File(in_dir + "p_i.xml") >> p_i
      self.model_inputs['p_i'] = p_i
      
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
      self.h_out << self.h
      self.pfo_out << self.pfo
      self.u_out << self.phi_solver.u
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


  # Write out xml checkpoint flies for h and phi
  def write_xmls(self, to_write = set([])):
    if len(to_write) == 0:
      File(self.check_dir + "h_" + str(self.check_index) + ".xml") << self.h
      File(self.check_dir + "phi_" + str(self.check_index) + ".xml") << self.phi
      File(self.check_dir + "u_" + str(self.check_index) + ".xml") << self.phi_solver.u
    else :
      if 'h' in to_write:
        File(self.check_dir + "h_" + str(self.check_index) + ".xml") << self.h
      if 'phi' in to_write:
        File(self.check_dir + "phi_" + str(self.check_index) + ".xml") << self.phi
      if 'u' in to_write:
        File(self.check_dir + "u_" + str(self.check_index) + ".xml") << self.phi_solver.u
        
    # Increment the checkpoint index
    self.check_index += 1
    
    
  # Updates the potentially time dependent expressions (m and u_b)
  def update_time_dependent_vars(self):
    # Check if the melt is an expression and if it has an attribute t. If so 
    # update t to the model time
    if isinstance(self.m, dolfin.Expression):
      if hasattr(self.m, 't'):
        # Update the time
        self.m.t = self.t  
        # Update expression form of m
        self.m_func.assign(project(self.m, self.V_cg))
        
    # Same deal for sliding speed
    if isinstance(self.u_b, dolfin.Expression):
      if hasattr(self.u_b, 't'):
        # Update the time
        self.u_b.t = self.t  
        # Update expression form of u_b
        self.u_b_func.assign(project(self.u_b, self.V_cg))