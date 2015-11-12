from dolfin import *
from dolfin_adjoint import *
from scipy.integrate import ode
from dolfin import MPI, mpi_comm_world
from colored import fg, attr

"""Solves ODE for the sheet height h with phi fixed."""

class HSolver():

  def __init__(self, model):
    
    # Process number
    self.MPI_rank = MPI.rank(mpi_comm_world())
    
    ### Get a few fields and parameters from the model
    
    # Effective pressure
    N = model.N   
    # Sliding speed
    u_b = model.u_b_func
    # Initial model time
    t0 = model.t
    # Rate factor
    A = model.pcs['A']
    # Distance between bumps
    l_r = model.pcs['l_r']
    # Bump height
    h_r = model.pcs['h_r']
    # Initial sheet height
    h0 = model.h.vector().array()
    

    ### Set up the sheet height ODE
      
    # Right hand side for the gap height ODE
    def rhs(t, h_n):
      # Ensure that the sheet height is positive
      h_n[h_n < 0.0] = 0.0
      # Sheet opening term
      w_n = u_b.vector().array() * (h_r - h_n) / l_r
      # Ensure that the opening term is non-negative
      w_n[w_n < 0.0] = 0.0
      # Sheet closure term
      v_n = A * h_n * N.vector().array()**3
      # Return the time rate of change of the sheet
      dhdt = w_n - v_n
      return dhdt
    
    # Set up ODE solver
    ode_solver = ode(rhs).set_integrator('vode',  method='adams', max_step = 60.0 * 5.0)
    ode_solver.set_initial_value(h0, t0)


    ### Set local variables    
    
    self.ode_solver = ode_solver
    self.model = model
    

  # Step the gap height h forward by dt
  def step(self, dt):
    # Use the current solution for h as the initial condition for the ODE solver
    # for this step. Normally, we could just use the solution from the previous
    # time step, but this is important if we want to reset the model (say manually
    # assign a different h)
    self.ode_solver.y[:] = self.model.h.vector().array()
    
    if self.MPI_rank == 0:
      print ('%sSolving for h... %s' % (fg(10), attr(0)))
      
    # Step h 
    self.ode_solver.integrate(self.model.t + dt)

    # Retrieve values from the ODE solver    
    self.model.h.vector().set_local(self.ode_solver.y)
    self.model.h.vector().apply("insert")
    
    if self.MPI_rank == 0:
      print ('%sDone. %s' % (fg(10), attr(0)))
  
 