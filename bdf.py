from dolfin import *
from dolfin import MPI, mpi_comm_world
import copy

class BDF(object):
  
  def __init__(self, F, U, U_dot, model):
    
    # Unknown at previous three time steps
    U1 = Function(model.V_cg)
    U2 = Function(model.V_cg)
    U3 = Function(model.V_cg)

    # Current time step
    h = Constant(1.0)
    # Last time step
    h1 = Constant(1.0)  
    # Time step before that
    h2 = Constant(1.0)
    
    # Time dependent weights for U + C1*U1 + C2*U2 = C3*F 
    w1 = h / h1
    C1 = -(1.0 + w1)**2 / (1.0 + 2.0*w1)
    C2 = w1**2 / (1.0 + 2.0*w1)
    C3 = h*((1.0 + w1)/(1.0 + 2.0*w1))
    
    
    ### Create first and second order forms
    
    # Approximation of derivative in first order form
    U_dot1 = (U - U1) / h
    # Approximation of derivative in second order form
    U_dot2 = (1.0/C3)*U + (C1/C3)*U1 + (C2/C3)*U2
    # First order form
    F1 = replace(F, {U_dot : U_dot1})
    # Second order form
    F2 = replace(F, {U_dot : U_dot2})
    # Forms for Jacobians    
    dU = TrialFunction(model.V_cg)
    J1 = derivative(F1, U, dU) 
    J2 = derivative(F2, U, dU) 
    
    
    ### Truncation error expression
    
    E1 = ((U - U1) / h)
    E2 = (1.0 + (h/h1))*((U1 - U2)/h1)
    E3 = (h/(h1*h2))*(U1 - U2)
    E = ((h + h1)/6.0)*(E1 - E2 + E3)   
    
    self.params = {}
    self.params['initial_step'] = 5.0 * 60.0    
    
    ### Set local vars
    self.U = U
    self.U1 = U1
    self.U2 = U2
    self.U3 = U3
    self.h = h
    self.h1 = h1
    self.F1 = F1
    self.F2 = F2
    self.J1 = J1
    self.J2 = J2
    self.model = model
    # Number of steps we've taken
    self.steps = 0
    # Function to store truncation error
    self.E = E
    self.E_func = Function(model.V_cg)
    
  
  def first_step(self, dt):
    # Step forward until we hit dt
    t = 0.0
    
    # Take two Euler steps
    h = self.params['initial_step']
    if h * 2.0 >= dt:
      h = dt / 3.0
    
    euler_steps = 0
    while euler_steps < 2:
      # Try to step forward with backward Euler and test convergence
      if self.try_step_d1(h):
        euler_steps += 1
        self.h1.assign(h)
        self.U1.assign(self.U)
        self.U2.assign(self.U1)
        t += h
      else :
        h /= 2.0
        
    # Now step forward with a second order method
    while t < dt:
      h = min(h, dt - t)
      
      if self.try_step_d2(h):
        
        # Estimate truncation error
        err = (self.U.array() - self.U2.array()) / float(self.h))
        err -= (1.0 + (float(self.h) / float(self.h1)))
        
        self.error.set_local()
        error = (float(h1) + h) / 6.0
        
        
        self.h1.assign(h)
        self.U1.assign(self.U)
        self.U2.assign(self.U1)
        t += h
        
        # Calculate new time step
        

      else :
        h /= 2.0
        
          
        
    
  # Try a step forward with backward Euler
  def try_step_d1(self, h):
    self.h.assign(h)
    solve(self.F1 == 0, self.U, self.model.d_bcs, J = self.J1, solver_parameters = self.model.newton_params)
    # Return true if it converges
    return True
    
  # Try a step forward second order BDF
  def try_step_d2(self, h):
    self.h.assign(h)
    solve(self.F2 == 0, self.U, self.model.d_bcs, J = self.J2, solver_parameters = self.model.newton_params)
    # Return true if it converges
    return True
  
    
    
    
    
    

    
