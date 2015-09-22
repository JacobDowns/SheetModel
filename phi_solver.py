from dolfin import *
from dolfin_adjoint import *

#parameters['form_compiler']['precision'] = 30

""" Solves phi with h fixed."""

class PhiSolver(object):
  
  def __init__(self, model):
    
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
    k = model.pcs['k']
    # Bump height
    h_r = model.pcs['h_r']
    # Distance between bumps
    l_r = model.pcs['l_r']
    # Exponents
    alpha = model.pcs['alpha']
    beta = model.pcs['beta']
    delta = model.pcs['delta']
    # Regularization parameter
    phi_reg = 1e-16
    

    ### Set up the PDE

    # Expression for effective pressure in terms of potential
    N = phi_0 - phi
    # Flux vector
    q = -Constant(k) * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(delta / 2.0) * grad(phi)
    # Opening term 
    w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)
    # Closing term
    v = Constant(A) * h * N**3
    
    # Function for storing the variation of the objective function 
    # Test function
    theta = TestFunction(model.V_cg)
    # Variational form for the PDE
    self.F = -dot(grad(theta), q) * dx + (w - v - m) * theta * dx
    # Get the Jacobian
    dphi = TrialFunction(model.V_cg)
    self.J = derivative(self.F, phi, dphi) 


    ### Solve the PDE once to create a record for dolfin-adjoint    
    phi.assign(phi_0)
    #self.solve_pde()
    
    
    ### Set up the variational principle

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
    J_phi = Constant((1. / beta) * k) * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(beta / 2.0) * dx  
    J_phi += (Constant(0.25 * A) * h * N**4) * dx
    J_phi += (w - m) * phi * dx
    
    self.J_hat = Functional(J_phi * dt[FINISH_TIME])
    self.rf = ReducedFunctional(self.J_hat, Control(phi, value = phi))
  

  # Solve the PDE for the potential
  def solve_pde(self):
    # Solve for potential
    solve(self.F == 0, self.model.phi, self.model.d_bcs, J = self.J, solver_parameters = self.model.newton_params)
    # Derive values form the potential
    self.model.update_phi()
    


  # Steps the potential forward
  def step(self):
    # Solve the pde
    #self.solve_pde()
    
    # If pressure is in the correct range then we're done. Otherwise we'll 
    # alter phi so it's in the correct range then use it as a starting guess
    # for the variational principle
    #over_or_under = self.phi_apply_bounds()
    
    #if over_or_under:
    minimize(self.rf, method = "L-BFGS-B", tol = 2e-08, bounds = (self.phi_min, self.phi_max), options = {"disp": True})
    
    quit()
    # Derive values from potential
    self.update_phi()
    
    plot(model.phi, interactive = True)
      
  
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
    indexes_over = phi_vals > phi_max_vals
    # Indexes that are underpressure
    indexes_under = phi_vals < phi_min_vals    
    
    if indexes_over.any() or indexes_under.any():
      # Correct under or overpressure values
      phi_vals[indexes_over] = phi_max_vals[indexes_over]
      phi_vals[indexes_under] = phi_min_vals[indexes_under]
  
      # Update phi    
      self.model.phi.vector().set_local(phi_vals)
      self.model.phi.vector().apply("insert")
      
      # If there were over or underpressure values return true
      return True
    
    # If pressure was in the correct range return false
    return False