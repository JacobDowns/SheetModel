from dolfin import *
from dolfin import MPI, mpi_comm_world
import numpy as np
import ufl
""" Steps forward (phi, h)"""
ufl.algorithms.apply_derivatives.CONDITIONAL_WORKAROUND = True

class Solver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # A reference to the model, which contains all the inputs we need
    self.model = model
    
    
    ### Get all the variables, inputs, and constants we need from the model
    
    # Combined (phi, h) unknown
    U = model.U
    # CG Function space
    V_cg = self.model.V_cg
    # Hydraulic potential 
    phi, h = split(U)
    # Time derivative of h
    h_prev = self.model.h_prev
    # Get melt rate
    m = model.m
    # Basal sliding speed
    u_b = model.u_b
    # Constraints
    phi_max = model.phi_max
    phi_min = model.phi_min
    # Potential at overburden pressure
    phi_0 = model.phi_0
    # Rate factor
    A = model.pcs['A']
    # Sheet conductivity
    k = model.k
    # Bump height
    h_r = model.h_r
    # Distance between bumps
    l_r = model.pcs['l_r']
    # Exponents
    alpha = model.pcs['alpha']
    delta = model.pcs['delta']
    # Regularization parameter
    phi_reg = 1e-16
    

    ### Set up the PDE for the potential 

    # Expression for effective pressure in terms of potential
    N = phi_0 - phi
    # Flux vector
    q = -k * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(delta / 2.0) * grad(phi)
    # Opening term 
    w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)
    #w = u_b * (h_r - h) / Constant(l_r)
    # Closing term
    v = Constant(A) * h * N**3
    # Test functions
    theta1, theta2 = TestFunctions(model.V)
    # Time step    
    dt = Constant(1.0)
    # PDE part of F
    F = -dot(grad(theta1), q)*dx + (w - v - m)*theta1*dx
    # ODE part of F
    F += ((h - h_prev)/dt - w + v)*theta2*dx
    # Get the Jacobian
    dU = TrialFunction(model.V)
    J = derivative(F, U, dU) 
    
                      
    # Set object variables                  
    self.F = F
    self.J = J
    self.U = U
    self.model = model
    self.q = q
    self.dt = dt
    
    
  # Step PDE for phi forward by dt. No constraints.
  def step(self, dt):
    self.dt.assign(dt)
    # Solve for potential
    solve(self.F == 0, self.U, self.model.d_bcs, J = self.J, solver_parameters = self.model.newton_params)
    
    self.model.update_U()
    self.h_prev.assign(self.model.h_cg)
    self.model.update_phi()      
