from dolfin import *
from dolfin import MPI, mpi_comm_world
import numpy as np

""" Solves phi with h fixed."""

class PhiSolver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # A reference to the model, which contains all the inputs we need
    self.model = model
    # Function space
    V_cg = model.V_cg
    self.V_cg  = V_cg
    # Get melt rate
    m = model.m
    # Sheet height
    h = model.h
    # Basal sliding speed
    u_b = model.u_b
    # Potential
    phi = model.phi
    self.phi = phi
    # Potential at previous time step
    phi1 = Function(V_cg)
    phi1.assign(phi)

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
    # Englacial void ratio
    e_v = model.pcs['e_v']
    # Density of water
    rho_w = model.pcs['rho_w']
    # Gravity
    g = model.pcs['g']
    # Exponents
    alpha = model.pcs['alpha']
    delta = model.pcs['delta']
    # Regularization parameter
    phi_reg = 1e-15

    

    ### Set up the PDE for the potential 

    # Expression for effective pressure in terms of potential
    N = phi_0 - phi
    # Flux vector
    q = -k * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(delta / 2.0) * grad(phi)
    # Opening term 
    w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)
    # Closing term
    v = Constant(A) * h * N**3
    # Constant in front of time derivative
    C = Constant(e_v/(rho_w * g))
    # Time step    
    dt = Constant(1.0)


    ### Variational form for backward Euler

    # Test function
    theta = TestFunction(model.V_cg)
    F = (theta*C*(phi - phi1) + dt*(-dot(grad(theta), q) + theta*(w - v - m)))*dx
    # Jacobian    
    dphi = TrialFunction(V_cg)
    J = derivative(F, phi, dphi) 
    
    prob = NonlinearVariationalProblem(F, phi, model.d_bcs, J)
    
    snes_solver_parameters = {"nonlinear_solver": "snes",
                      "snes_solver": {"linear_solver": "lu",
                                      "maximum_iterations": 30,
                                      "report": True,
                                      "line_search" : 'basic',
                                      "error_on_nonconvergence": False, 
                                      "relative_tolerance" : 1e-10,
                                      "absolute_tolerance" : 1e-6}}
    
    solver = NonlinearVariationalSolver(prob)
    solver.parameters.update(snes_solver_parameters)
    
    self.phi = phi  
    self.phi1 = phi1
    self.dt = dt
    self.F = F
    self.J = J
    self.model = model
    self.q = q
    self.solver = solver
 
  
  # Step PDE for phi forward by dt. No constraints.
  def step(self, dt):
    # Assign time step
    self.dt.assign(dt)
    # Solve for potential
    #(i, converged) = self.solver.solve()
    solve(self.F == 0, self.phi, self.model.d_bcs, J = self.J, solver_parameters = self.model.newton_params)
    # Update phi1
    self.phi1.assign(self.phi)
    # Update fields derived from phi
    self.model.update_phi()
    
