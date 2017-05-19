from dolfin import *
from dolfin import MPI, mpi_comm_world

""" Steps forward (phi, h)"""

class Solver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # A reference to the model, which contains all the inputs we need
    self.model = model
    # Combined (phi, h) unknown
    U = U
    # CG Function space
    V_cg = self.model.V_cg
    # Hydraulic potential and sheet height 
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
    # Closing term
    v = Constant(A) * h * N**3
    
    # Test function
    theta1, theta2 = TestFunction(model.V_cg)
    # PDE part of F
    F = -dot(grad(theta1), q)*dx + (w - v - m)*theta1*dx
    # Time step    
    dt = Constant(1.0)
    # ODE part of F
    F += ((h - h_prev)/dt - w + v)*theta2*dP
    # Get the Jacobian
    dphi = TrialFunction(model.V_cg)
    J = derivative(F, phi, dphi) 

                      
    # Set object variables                  
    self.F = F
    self.J = J
    self.phi = phi
    self.model = model
    self.q = q
    
    
  # Step PDE for phi forward by dt. No constraints.
  def step(self):
    # Solve for potential
    solve(self.F == 0, self.phi, self.model.d_bcs, J = self.J, solver_parameters = self.model.newton_params)
    # Update phi
    self.model.update_phi()  
    
    
  # Step PDE for phi forward by dt. Constrain using SNES solver. 
  def step_constrained(self):
    # Solve for potential
    (i, converged) = self.phi_solver.solve()
    # Update phi
    self.model.update_phi()  
    
