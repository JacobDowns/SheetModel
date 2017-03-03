from dolfin import *
from dolfin import MPI, mpi_comm_world
import numpy as np
from bdf_helper import *

class BDF(object):
  
  def __init__(self, F, y, y_prime, bcs = [], params = None):
    
    # Process number
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # Get the function space
    V  = y.function_space()
    
    #Form F(Y_prime, Y, t) = 0
    self.F = F
    
    # Unknown at previous three time steps
    # Form is [y0, ..., y5] with y0 most recent
    ys = [y]
    for i in range(1,6):
      ys.append(Function(V, name = 'y_' + str(i)))
      
    # Set an initial condition for degree 1 method (backward Euler)
    ys[1].assign(y)
    
    # Times of last 6 solutions
    # Form is [t0, ..., t5] with t0 most recent
    ts = np.zeros(6)
    
  
    ### Create form for computing Jacobian of F used in Netwon's method
  
    """
    Linearizing F for Netwon's method we want Fenics to generate the matrix
      
      F_y + (alpha_0 / h) * F_(y')
      
    Here F_y is the derivative of F wrt Y, F_(y') is derivative wrt y', alpha_0
    is a constant that depends on the order and h is the time step. Here we'll
    create a form that generates this matrix. 
    """
    
    # Leading coefficients for each order
    alpha_0s = np.array([1.0, 3.0/2.0, 11.0/6.0, 25.0/12.0, 137.0/60.0])
    # Standard trial function
    dy = TrialFunction(V)
    # Shift is a stand in for (alpha_0 / h). It depends on the time step and
    # the method order
    shift = Constant(1.0)
    # F_y
    F_y_form = derivative(F, y, dy)
    # F_(y')
    F_y_prime_form = derivative(F, y_prime, dy)
    # Jocobian 
    J = F_y_form + shift*F_y_prime_form
    
    
    ### Predicted values
    
    self.y_p = Function(V)
    self.y_prime_p = Function(V)
      
    
    ### BDF solver params
    
    self.params = params
    if self.params == None:
      self.params = {}
      self.params['initial_step'] = 0.0005 
      self.params['max_tries'] = 3
      self.params['tol'] = 1e-8
      self.params['lte_tries'] = 10
      self.params['newton_atol'] = 1e-7
      self.params['newton_rtol'] = 1e-11
      self.params['newton_max_iters'] = 15
      
    
    ### Create some attributes
    
    self.ys = ys
    self.ts = ts
    # Create a helper object to help compute coefficients and LTE 
    self.bdf_helper = BDF_Helper()
    # Set to true once the method has been bootstrapped
    self.bootstrapped = False
    # This function will store LTE error 
    self.lte_func = Function(V)
      
    
  # Take a step with BDFk
  def step_dk(self, k, h):

    success = False
    tries = 0 
    max_tries = self.params['max_tries']
    
    # Get the BDFk solver
    solver = self.solvers[k]
    # In order to compute coefficients for degree k, bdf_helper needs the times
    # of previous solutions as well as the time we want to advance to. Here
    # we're proposing to advance to time t0 + h with previous times t0, t1, ...
    times = np.zeros(6)
    times[1:] = self.ts[0:5]
    
    # Try different time steps until the solver converges. If the solver doesn't 
    # converge after halving the time step some number of times, then just accept the solution.
    while not success and tries < max_tries:
      times[0] = self.ts[0] + h
      
      # Use our friendly bdf helper to compute some coefficients!
      self.bdf_helper.ts = times
      cofs_k = self.bdf_helper.compute_coefficients(k)[:(k+1)]

      # Now we need to update the coefficients in the BDFk form
      for i in range(k+1):
        self.cofs[k][i].assign(cofs_k[i])
      
      # Now that the coefficients have been set, solve
      (i, converged) = solver.solve()
      success = converged 
      tries += 1

      if success:
        break
      else:
        # If solver didn't converge cut time step in half
        h /= 2.0
    
    if tries == max_tries and self.MPI_rank == 0:
      print "Solver did not converge after " + str(tries) + " tries. Accepting solution as is."
      
    # Return the time step we took
    return h, times
    
  
  # Try steps with BDFk of different sizes until the LTE is below the threshold
  def step_lte_dk(self, k, h):    
    for i in range(self.params['lte_tries']):
      # Try a step h 
      h, times = self.step_dk(k, h)
      # Check LTE 
      lte_norm = self.compute_lte(k, times)
      # Error tolerance
      etol = self.params['tol']
      print etol
      # Propose a new time step based on error
      alpha = (((0.9) * etol) / lte_norm)**(1.0 / (k + 1.0))
      hnew = alpha * h
      
      print (0.9 * etol)
      print lte_norm
      print alpha
      print "lte_norm", lte_norm
      print "hnew", hnew
      
      quit()
      
      # Check if we should accept the step
      if lte_norm * dt <= self.params['tol']:
        self.accept_step(h)
        return hnew
      else :
        # If error is too much try a smaller time step
        h = hnew
        
    print "Couldn't reduce LTE error sufficiently in " + self.params['lte_tries'] + " tries. Accepting solution as is."
    quit()
  
  
  # Initialize the ODE solver
  def bootstrap(self, h):
    
    for k in range(1,4):
      
      # Initially take two steps with backward Euler so we can compute LTE
      if k == 1:
        h, times = self.step_dk(k, h)
        self.accept_step(h)
        
        self.step_lte_dk(k, h)
        
      else : 
        # Take a step with a method of degree k
        pass
        
      quit()
      # Compute a local truncation error
      self.compute_lte(k, times)
        

  # Accept a time step from t0 to t0 + h     
  def accept_step(self, h):
    # Update times
    self.ts = np.roll(self.ts, 1)
    self.ts[0] = self.ts[1] + h
    # Update solutions 
    for i in range(5,0,-1):
      self.Ys[i].assign(self.Ys[i-1])
      

  # Compute the local truncation error of order form BDFk. LTE is a linear combination
  # of Y0, ... , Y5
  def compute_lte(self, k, times):
    print times
    self.bdf_helper.ts = times
    # Compute coefficients
    lte_cofs = self.bdf_helper.compute_lte_coefficients(k)
    # Compute LTE
    lte = lte_cofs[0] * self.Ys[0].vector().array()
    for i in range(1, k + 2):
      lte += lte_cofs[i] * self.Ys[i].vector().array()

    # Pute LTE in a function      
    self.lte_func.vector().set_local(lte)
    self.lte_func.vector().apply("insert")
    
    # compute norm
    return self.lte_func.vector().norm('linf')
    
    
  
  def step_d1(self, h):
    # Take a step
    self.try_step(h, self.d1_solver) 
    
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
    
    
  # Based on past values of the solution we can try to predict what the solution
  # and its derivative will be at t+h
    

  # Special Newton solver 
  def newton_solve(self, h):
    # We assume here that shift has been correctly for the method order
    
    i = 0
    
    # Convergence params
    atol = self.params['newton_atol']
    rtol = self.params['newton_rtol']
    max_iters = self.prams['newton_max_iters']
    a_err = 2.0*atol    
    r_err = 2.0*rtol
    
    while i < max_iters and a_err > atol and r_err > rtol :
      pass
      
  
    
    
    
    
    

    
