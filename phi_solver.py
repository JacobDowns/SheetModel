from dolfin import *
from dolfin_adjoint import *

""" Solves phi with h fixed."""

class PhiSolver(object):
  
  def __init__(self, model):
    

    # Get melt rate
    m = model.m
    # Sheet height
    h = model.h
    # Basal sliding speed
    u_b = model.u_b
    # Potential
    phi = model.phi
    # Starting guess for phi
    phi_init = model.phi_init
    # Potential at 0 pressure
    phi_m = model.phi_m
    # Potential at overburden pressure
    phi_0 = model.phi_0
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
    # Regularization parameter
    phi_reg = 1e-15
    

    ### Set up the potential functional

    # Expression for effective pressure in terms of potential
    N = phi_0 - phi
    # Opening term 
    w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)

    # Define some upper and lower bounds for phi
    phi_min = Function(model.V_cg)
    phi_max = Function(model.V_cg)
    
    phi_min.assign(phi_m)
    phi_max.assign(phi_0)
    
    # On Dirichlet boundary conditions set the min and the max to be equal
    for bc in model.d_bcs:
      bc.apply(phi_min.vector())
      bc.apply(phi_max.vector())
    
    # Make sure that the initial guess for phi is within the constraints
    phi_vals = phi_init.vector().array()
    phi_max_vals = phi_max.vector().array()
    phi_min_vals = phi_min.vector().array()
    
    indexes_over = phi_vals > phi_max_vals
    indexes_under = phi_vals < phi_min_vals
    phi_vals[indexes_over] = phi_max_vals[indexes_over]
    phi_vals[indexes_under] = phi_min_vals[indexes_under]
    
    phi_init.vector().set_local(phi_vals)
    phi_init.vector().apply("insert")
    #phi.assign(phi_init)
    
    
    # Functional for the potential
    J_phi = Constant((1. / beta) * k) * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(beta / 2.0) * dx  
    J_phi += (Constant(0.25 * A) * h * N**4) * dx
    J_phi += (w - m) * phi * dx
    
    # Add any non-zero Neumann boundary conditions. m is the boundary marker
    # number for the boundary facet function. c is the flux value. 
    #for (m, c) in model.n_bcs:
    #  J_phi += Constant(c) * theta * ds(m)
    
    J_hat = Functional(J_phi * dt[FINISH_TIME])
    reduced_functional = ReducedFunctional(J_hat, Control(phi))
    
    
    ### Assign local variables

    self.reduced_functional = reduced_functional
    self.model = model
    self.phi_min = phi_min
    self.phi_max = phi_max
    self.J_hat = J_hat


  # Steps the potential forward
  def step(self):
    minimize(self.reduced_functional, method = "L-BFGS-B", tol = 2e-08, bounds = (self.phi_min, self.phi_max), options = {"disp": True})
    
    # Derive values from potential
    self.update_phi()
      
  
  