from dolfin import *
from dolfin_adjoint import *
from dolfin import MPI, mpi_comm_world
from colored import fg, attr

parameters['form_compiler']['precision'] = 30

""" Solves phi with h fixed."""

class PhiSolver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
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
    # Potential at 0 pressure
    phi_m = model.phi_m
    # Potential at overburden pressure
    phi_0 = model.phi_0
    # This is necessary due to a bizarre dolfin-adjoint bug. For some reason
    # reprojecting gets rid of it. 
    phi_0 = project(phi_0, model.V_cg)
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
    u = Function(model.V_cg)
    u.assign(phi_m)
    self.u = u
    # Expression for effective pressure in terms of potential
    Nu = phi_0 - u
    # Flux vector
    qu = -k * h**alpha * (dot(grad(u), grad(u)) + phi_reg)**(delta / 2.0) * grad(u)
    # Opening term 
    wu = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)
    # Closing term
    vu = Constant(A) * h * Nu**3
    
    # Function for storing the variation of the objective function 
    # Test function
    theta = TestFunction(model.V_cg)
    # Variational form for the PDE
    self.F = -dot(grad(theta), qu) * dx + (wu - vu - m) * theta * dx
    # Get the Jacobian
    du = TrialFunction(model.V_cg)
    self.J = derivative(self.F, u, du) 

  
    ### Set up the variational principle
    
    phi.assign(phi_0)
    # Expression for effective pressure
    N = phi_0 - phi
    # Opening term 
    w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)

    # Define some upper and lower bounds for phi
    self.phi_min = Function(model.V_cg)
    self.phi_max = Function(model.V_cg)
    self.phi_min.assign(phi_m)
    self.phi_max.assign(phi_0)
    
    # On Dirichlet boundary conditions set the min and the max to be equal
    for bc in model.d_bcs:
      bc.apply(self.phi_min.vector())
      bc.apply(self.phi_max.vector())    
    
    # Functional for the potential
    J_phi = Constant(1. / beta) * k * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(beta / 2.0) * dx  
    J_phi += (Constant(0.25 * A) * h * N**4) * dx
    J_phi += (w - m) * phi * dx
    
    self.J_hat = Functional(J_phi * dt[FINISH_TIME])
    self.rf = ReducedFunctional(self.J_hat, Control(phi, value = phi))
  

  # Internal function that solves the PDE for u
  def __solve_pde__(self):
    # Solve for potential
    solve(self.F == 0, self.u, self.model.d_bcs, J = self.J, solver_parameters = self.model.newton_params)    

    
  # External function that solves PDE for u then copies it to model.phi and updates
  # any model fields related to phi
  def solve_pde(self):
    self.__solve_pde__()
    # Copy PDE solution u to model phi
    self.model.phi.assign(self.u)    
    # Update phi
    self.model.update_phi()


  # Internal function that solves optimization problem for model.phi. Tolerance
  # is the usual solver tolerance and scale rescales the problem. Increasing the
  # scale can prevent premature convergence
  def __solve_opt__(self, tol, scale):
    # Solve the optimization problem
    minimize(self.rf, method = "L-BFGS-B", scale = scale, tol = tol, bounds = (self.phi_min, self.phi_max), options = {"disp": True})


  # External function that solves optimization problem for model.phi then updates 
  # any fields related to phi 
  def solve_opt(self):
    scale = self.model.opt_params['scale']
    tol = self.model.opt_params['tol']
    
    self.__solve_opt__(tol, scale)
    self.model.update_phi()
    
    parameters["adjoint"]["stop_annotating"] = True
    

  # Steps the potential forward with h fixed
  def step(self):
    # Pretty colors! (:
    #print('%s' % (fg(231)))
    
    # Solve the PDE with an initial guess of phi_m
    self.u.assign(self.model.phi_m)
    self.__solve_pde__()
    
    # Copy the solution u to phi
    self.model.phi.assign(self.u)    
    self.model.update_phi()
    
    # Check if there is any over or under pressure on this process
    local_over_or_under = self.phi_apply_bounds()
    # This will be 1 if there is over or underpressure on any process and 0
    # otherwise
    global_over_or_under = MPI.max(mpi_comm_world(), local_over_or_under)        
    
    # If we do get over or under pressure, we'll solve the optimization problem
    # to correct it
    if global_over_or_under:
      scale = self.model.opt_params['scale']
      tol = self.model.opt_params['tol']
      # Solve the optimization problem with the PDE solution as an initial guess
      self.__solve_opt__(tol, scale)
    
    # Update any fields derived from phi
    self.model.update_phi()
    
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