from dolfin import *
from scipy.integrate import ode
import numpy as np

"""Solves for h and S with phi fixed."""

class HSSolver():

  def __init__(self, model):
    
    ### Get a few fields and parameters from the model
    
    self.model = model
    # Effective pressure
    N = model.N   
    # Sheet height on edges
    h_cr = model.h_cr
    # Effective pressure on edges
    N_cr = model.N_cr
    # Derivative of phi over edges    
    dphi_ds_cr = model.dphi_ds_cr
    # Sliding speed
    u_b = model.u_b
    # Sheet conductivity
    k_cr = model.k_cr
    # Bump height
    h_r = model.h_r
    # Initial model time
    t0 = model.t
    # Rate factor
    A = model.pcs['A']
    # Distance between bumps
    l_r = model.pcs['l_r']
    # Density of ice
    rho_i = model.pcs['rho_i']
    # Latent heat
    L = model.pcs['L']
    # Channel conductivity
    k_c = model.pcs['k_c']
    # Sheet width under channel
    l_c = model.pcs['l_c']
    # Exponent
    alpha = model.pcs['alpha']
    delta = model.pcs['delta']
    # Regularization parameter
    phi_reg = 1e-16   
    
    
    ### Static arrays used in the ODE rhs
    
    # Mask used to prevent opening on exterior edges
    local_mask = model.mask.vector().array()
    # Initial sheet height
    h0 = model.h.vector().array()
    # Initial channel areas
    S0 = model.S.vector().array()
    # Length of h vector
    h_len = len(h0) 
    # Bump height vector
    h_r_n = h_r.vector().array()
    

    ### Set up the sheet height and channel area ODEs
    
    # Right hand side for the gap height ODE
    def h_rhs(t, h_n) :
      # Ensure that the sheet height is positive
      h_n[h_n < 0.0] = 0.0
      # Sheet opening term
      w_n = u_b.vector().array() * (h_r_n - h_n) / l_r
      # Ensure that the opening term is non-negative
      w_n[w_n < 0.0] = 0.0
      # Sheet closure term
      v_n = A * h_n * N.vector().array()**3
      # Return the time rate of change of the sheet
      dhdt = w_n - v_n
      return dhdt
      
    # Right hand side for the channel area ODE
    def S_rhs(t, S_n):
      # Ensure that the channel area is positive
      S_n[S_n < 0.0] = 0.0
      # Get effective pressures, sheet thickness on edges.
      N_n = N_cr.vector().array()
      # Get midpoint values of sheet thickness
      h_n = h_cr.vector().array()
      # Array form of the derivative of the potential 
      phi_s = dphi_ds_cr.vector().array()  
      # Along channel flux
      Q_n = -k_c * S_n**alpha * abs(phi_s + phi_reg)**delta * phi_s
      # Flux of sheet under channel
      q_n = k_cr.vector().array() * h_n**alpha * abs(phi_s + phi_reg)**delta * phi_s
      # Dissipation melting due to turbulent flux
      Xi_n = abs(Q_n * phi_s) + abs(l_c * q_n * phi_s)
      # Creep closure
      v_c_n = A * S_n * N_n**3
      # Total opening rate
      v_o_n = Xi_n / (rho_i * L)
      # Dissalow negative opening rate where the channel area is 0
      #v_o_n[v_o_n[S_n == 0.0] < 0.0] = 0.0
      # Calculate rate of channel size change
      dsdt = local_mask * (v_o_n - v_c_n)
      return dsdt
      
    # Combined right hand side for h and S
    def rhs(t, Y):
      Ys = np.split(Y, [h_len])
      h_n = Ys[0]
      S_n = Ys[1]
      
      dhdt = h_rhs(t, h_n)
      dsdt = S_rhs(t, S_n)
      
      return np.hstack((dhdt, dsdt))
    
    # ODE solver initial condition
    Y0 = np.hstack((h0, S0))
    # Set up ODE solver
    ode_solver = ode(rhs).set_integrator('vode', method = 'adams', max_step = 60.0 * 5.0)
    ode_solver.set_initial_value(Y0, t0)


    ### Set local variables    
    
    self.ode_solver = ode_solver
    self.model = model
    self.h_len = h_len
    

  # Step  h and S forward by dt
  def step(self, dt):
    # Step h and S forward
    self.ode_solver.integrate(self.model.t + dt)

    # Retrieve values from the ODE solver    
    Y = np.split(self.ode_solver.y, [self.h_len])
    self.model.h.vector().set_local(Y[0])
    self.model.h.vector().apply("insert")
    self.model.S.vector().set_local(Y[1])
    self.model.S.vector().apply("insert")
    
    # Update any fields derived from h or S
    self.model.update_S()
    self.model.update_h()
    
  
  
 