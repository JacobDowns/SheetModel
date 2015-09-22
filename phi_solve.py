from dolfin import *
from dolfin_adjoint import *
from constants import *
from dolfin import MPI, mpi_comm_world

# Variational potential solver
class VPS(object):
  
  def __init__(self, model):
    
    # CG function space
    V_cg = FunctionSpace(model.mesh, "CG", 1)    
  
    ### Sheet model 
    
    # Effective pressure
    #N = phi_0 - phi
    # Sheet opening term
    #w = u_b * (Constant(h_r) - h) / Constant(l_r)
    # Initial guess for phi
    #phi.assign(phi_0)
    
    u = Function(V_cg)


    ### Functional for the potential ###
    
    #J1 = Constant((1.0 / beta) * k) * h**alpha *(dot(grad(phi), grad(phi)) + phi_reg)**(beta / 2.0)
    #J2 = Constant(0.25 * A) * h * N**4 
    #J3 = (w - m) * phi 
    
    # The full functional as a fenics form 
    J_phi = u * dx
    J = Functional(J_phi * dt[FINISH_TIME])
    J_hat = ReducedFunctional(J, Control(u, value = u))
    
    # Upper and lower bounds for the potential    
    #phi_min.assign(phi_m)
    #phi_max.assign(phi_0)
    #bc.apply(phi_max.vector())
    
    
    minimize(J_hat, method = "L-BFGS-B", tol = 2e-08, options = {"disp": True})
    

  def solve(self):    
    m_opt = minimize(self.J_hat, method = "TNC", tol = 1.05e-7, bounds = (self.phi_min, self.phi_max), options = {"disp": True})
      
     

  