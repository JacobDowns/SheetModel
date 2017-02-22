from dolfin import *
from dolfin import MPI, mpi_comm_world
import copy

class BDF(object):
  
  def __init__(self, F, U, U_dot, model):
    
    # Process number
    self.MPI_rank = MPI.rank(mpi_comm_world())
    
    # Unknown at previous three time steps
    U1 = Function(model.V_cg)
    U1.assign(U)
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
    
    problem_d1 = NonlinearVariationalProblem(F1, U, model.d_bcs, J1)
    problem_d2 = NonlinearVariationalProblem(F2, U, model.d_bcs, J2)
    
    snes_solver_parameters = {"nonlinear_solver": "snes",
                      "snes_solver": {"linear_solver": "lu",
                                      "maximum_iterations": 30,
                                      "report": True,
                                      "line_search" : 'basic',
                                      "error_on_nonconvergence": False, 
                                      "relative_tolerance" : 1e-10,
                                      "absolute_tolerance" : 1e-6}}
    
    d1_solver = NonlinearVariationalSolver(problem_d1)
    d1_solver.parameters.update(snes_solver_parameters)
    d2_solver = NonlinearVariationalSolver(problem_d2)
    d2_solver.parameters.update(snes_solver_parameters)
    

    
    
    ### Truncation error expression
    
    #E1 = ((U - U1) / h)
    #E2 = (1.0 + (h/h1))*((U1 - U2)/h1)
    #E3 = (h/(h1*h2))*(U1 - U2)
    #E = ((h + h1)/6.0)*(E1 - E2 + E3)   
    
    self.params = {}
    self.params['initial_step'] = 10.0 * 60.0   
    self.params['max_tries'] = 3
    
    ### Set local vars
    self.U = U
    self.U1 = U1
    self.U2 = U2
    self.U3 = U3
    self.h = h
    self.h1 = h1
    self.h2 = h2
    self.d1_solver = d1_solver
    self.d2_solver = d2_solver
    self.model = model
    self.w1 = w1
    # Number of steps we've taken
    self.steps = 0
    # Function to store truncation error
    #self.E = E
    #self.E_func = Function(model.V_cg)
  
    
    self.bootstrapped = False
   
  
  def step(self, dt):
    t = 0.0
    
    # If the ODE solver isn't primed then do it now
    if not self.bootstrapped:
      # Make sure taking two priming steps doesn't take us beyond dt
      h_init = self.params['initial_step']
      if h_init * 2.0 >= dt:
        h_init = dt / 3.0
      
      # Take degree 1 step
      self.step_d1(h_init)
      t += float(self.h)
      
      # Take degree 2 step
      self.try_step(h_init, self.d2_solver) 
      t += float(self.h)
      
      # Update fields
      self.h2.assign(self.h1)
      self.h1.assign(self.h)
      self.U3.assign(self.U2)
      self.U2.assign(self.U1)
      self.U1.assign(self.U)
    
    # Take steps until we reach dt
    while t < dt:
      h = min(h, dt - t)
      
     
        
          
    
    
  # Take a step with the degree 1 or 2 BDF
  def try_step(self, h, problem):
    # Try different time steps until the solver converges. If the solver doesn't 
    # converge after halving the time step some number of times, then just accept the solution.
    success = False
    tries = 0 
    max_tries = self.params['max_tries']
    while not success and tries < max_tries:
      #print "tries", tries
      #print "h"
      self.h.assign(h)
      
      #print "w1", float(self.w1)
      (i, converged) = problem.solve()
      success = converged 
      tries += 1

      if success:
        break
      else:
        # If solver didn't converge cut time step in half
        h /= 2.0
    
    if tries == max_tries and self.MPI_rank == 0:
      print "Solver did not converge after " + str(tries) + " tries. Accepting solution as is."
  
  
  def step_d1(self, h):
    
    # Take a step
    #self.try_step(h, self.d1_solver)
    (i, converged) = self.d1_solver.solve()   
    
    File('U.pvd') << self.U
    
    # Update fields
    self.h2.assign(self.h1)
    self.h1.assign(self.h)
    self.U3.assign(self.U2)
    self.U2.assign(self.U1)
    self.U1.assign(self.U)
    
    
  # Try a step forward second order BDF
  def try_step_d2(self, h):
    self.h.assign(h)
    solve(self.F2 == 0, self.U, self.model.d_bcs, J = self.J2, solver_parameters = self.model.newton_params)
    # Return true if it converges
    return True
  
    
    
    
    
    

    
