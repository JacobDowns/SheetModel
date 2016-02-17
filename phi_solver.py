from dolfin import *
from dolfin import MPI, mpi_comm_world
from colored import fg, attr
from scipy.optimize import fmin_l_bfgs_b
import numpy as np

parameters['form_compiler']['precision'] = 30

""" Solves phi with h fixed."""

class PhiSolver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # Local to global dof map    
    self.global_indexes = model.V_cg.dofmap().dofs()
    # Array of all indexes. This is used to gather values on multiple processes
    # into a single vector
    self.indexes = np.array(range(model.V_cg.dim()), dtype = np.intc)
    # A reference to the model, which contains all the inputs we need
    self.model = model
    # Get melt rate
    m = model.m
    # Sheet height
    h = model.h
    # Basal sliding speed
    u_b = model.u_b
    # Potential
    phi = model.phi
    self.phi = phi
    # Potential at 0 pressure
    phi_m = model.phi_m
    # Potential at overburden pressure
    phi_0 = model.phi_0
    # Rate factor
    A = model.pcs['A']
    # Sheet conductivity
    k = model.k
    # Bump height
    h_r = model.pcs['h_r']
    # Distance between bumps
    l_r = model.pcs['l_r']
    # Exponents
    alpha = model.pcs['alpha']
    beta = model.pcs['beta']
    delta = model.pcs['delta']
    # Regularization parameter
    phi_reg = 1e-15
    

    ### Set up the PDE for the potential 

    # Unknown 
    phi.assign(phi_m)
    # Expression for effective pressure in terms of potential
    N = phi_0 - phi
    # Flux vector
    q = -k * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(delta / 2.0) * grad(phi)
    # Opening term 
    w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)
    # Closing term
    v = Constant(A) * h * N**3
    
    # Test function
    theta = TestFunction(model.V_cg)
    # Variational form for the PDE
    self.F = -dot(grad(theta), q) * dx + (w - v - m) * theta * dx
    # Get the Jacobian
    dphi = TrialFunction(model.V_cg)
    self.J = derivative(self.F, phi, dphi) 

  
    ### Set up the variational principle for determining potential 
  
    # Functional for the potential
    self.J_phi = Constant(1. / beta) * k * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(beta / 2.0) * dx  
    self.J_phi += (Constant(0.25 * A) * h * N**4) * dx
    self.J_phi += (w - m) * phi * dx
    
    
    # Function to store the variation of the functional
    self.dJ_phi = Function(model.V_cg)
    # Vector to store the global variation of the functional 
    self.dJ_phi_global = Vector()
    # Vector to store global value of phi
    self.phi_global = Vector()

    # Define some upper and lower bounds for phi

    # Lower bound as function
    self.phi_min = Function(model.V_cg)
    self.phi_min.assign(phi_m)
    # Upper bound as function    
    self.phi_max = Function(model.V_cg)
    self.phi_max.assign(phi_0)
    
    # On Dirichlet boundary conditions set the min and the max to be equal
    for bc in model.d_bcs:
      bc.apply(self.phi_min.vector())
      bc.apply(self.phi_max.vector())  

    # Lower bound as a vector
    self.phi_min_global = Vector()
    phi_m.vector().gather(self.phi_min_global, self.indexes)
    # Upper bound as a vector
    self.phi_max_global = Vector()
    phi_0.vector().gather(self.phi_max_global, self.indexes)
    
    # Bounds for lbfgsb
    self.bounds = zip(self.phi_min_global.array(), self.phi_max_global.array())
    
    
  # Python function version of the functional
  def __J_phi_func__(self, x):
    # Assign phi
    self.__set_phi__(x)
    
    #plot(self.phi, interactive = True)
    
    # Function value of functional
    F =  assemble(self.J_phi)
    
    # Scale and return the value
    return self.model.opt_params['scale'] * F
    
    
  # Python function for the variation of the functional
  def __var_J_phi_func__(self, x):
    # Assign phi
    self.__set_phi__(x)
    # Store the variation of the functional in a function
    self.dJ_phi.vector().set_local(assemble(self.F).array())
    self.dJ_phi.vector().apply("insert")
    
    # Enforce 0 variation on Dirichlet boundaries 
    #for bc in self.model.d_bcs:
    #  bc.apply(self.dJ_phi.vector())
    
    # Now gather the variation into a single global vector and return it
    self.dJ_phi.vector().gather(self.dJ_phi_global, self.indexes)
    # Save the array     
    dF = np.array(self.dJ_phi_global.array())
    # Reset the vector or petsc complains
    self.dJ_phi_global = Vector()
    
    return self.model.opt_params['scale'] * dF


  # Sets the local part of phi given a global vector x
  def __set_phi__(self, x):
    self.model.phi.vector().set_local(x[self.global_indexes])
    self.model.phi.vector().apply("insert")


  # Internal function that solves the PDE for u
  def __solve_pde__(self):
    # Solve for potential
    solve(self.F == 0, self.phi, self.model.d_bcs, J = self.J, solver_parameters = self.model.newton_params)    

    
  # External function that solves PDE for u then copies it to model.phi and updates
  # any model fields related to phi
  def solve_pde(self):
    self.__solve_pde__() 
    # Update phi
    self.model.update_phi()
    

  # Internal function that solves optimization problem for model.phi. Tolerance
  # is the usual solver tolerance and scale rescales the problem. Increasing the
  # scale can prevent premature convergence
  def __solve_opt__(self):
    self.phi.assign(self.phi_max)
    # Tolerance
    tol = self.model.opt_params['tol']
    # Gather phi into a single vector called phi_global
    self.phi.vector().gather(self.phi_global, self.indexes)
    # Display first processor only
    disp = (self.MPI_rank == 0)
    
    # Initial condition
    phi_start = np.array(self.phi_global.array())
    # Solve the optimization problem
    phi_opt, f, d = fmin_l_bfgs_b(self.__J_phi_func__, phi_start, self.__var_J_phi_func__, 
                  pgtol = tol, disp = disp)
                  
    # Now assigm phi_opt to phi
    self.__set_phi__(phi_opt)

    plot(self.phi, interactive = True)


  # External function that solves optimization problem for model.phi then updates 
  # any fields related to phi 
  def solve_opt(self):    
    self.__solve_opt__(tol, scale)
    self.model.update_phi()
    

  # Steps the potential forward with h fixed
  def step(self):
    # Pretty colors! (:
    #print('%s' % (fg(231)))
    
    # Solve the PDE with an initial guess of phi_m
    self.phi.assign(self.model.phi_m)
    self.__solve_pde__()
    # Update fields related to phi
    self.model.update_phi()
    
    
    # Check if there is any over or under pressure on this process
    local_over_or_under = self.phi_apply_bounds()
    # This will be 1 if there is over or underpressure on any process and 0
    # otherwise
    global_over_or_under = MPI.max(mpi_comm_world(), local_over_or_under)        
    
    # If we do get over or under pressure, we'll solve the optimization problem
    # to correct it
    if global_over_or_under:
      # Solve the optimization problem with the PDE solution as an initial guess
      self.__solve_opt__()
    
    # Update any fields derived from phi
    self.model.update_phi()
    
    plot(self.model.pfo, interactive = True)
    
    # No more pretty colors. ):
    #print('%s' % (attr(0))),
      
  
  # Correct the potential so that it is above 0 pressure and below overburden.
  # Return True if underpressure was present?
  def phi_apply_bounds(self):
    
    # Array of values for phi
    phi_vals = self.model.phi.vector().array()
    # Array of minimum values
    phi_max_vals = self.phi_max.vector().array()
    # Array of maximum values
    phi_min_vals = self.phi_min.vector().array()
    
    # Indexes in the array of phi vals that are overpressure
    indexes_over = phi_vals > phi_max_vals + 1e-3
    # Indexes that are underpressure
    indexes_under = phi_vals < phi_min_vals - 1e-3    
    
    phi_vals[indexes_over] = phi_max_vals[indexes_over]
    phi_vals[indexes_under] = phi_min_vals[indexes_under]
  
    # Update phi    
    self.model.phi.vector().set_local(phi_vals)
    self.model.phi.vector().apply("insert")
    
    # If there were over or underpressure values return true
    if indexes_over.any() or indexes_under.any():
      return True
    
    # If pressure was in the correct range return false
    return False