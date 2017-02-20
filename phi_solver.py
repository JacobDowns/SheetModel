from dolfin import *
from dolfin import MPI, mpi_comm_world
from petsc4py import PETSc
import numpy as np
from bdf import *

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
    # Time derivative of phi
    phi_dot = Function(V_cg)
    self.phi_dot = phi_dot
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
    # Test function
    theta = TestFunction(model.V_cg)
    # Variational form for the PDE
    F = -dot(grad(theta), q) * dx + (w - v - m) * theta * dx
    # Get the Jacobian
    dphi = TrialFunction(model.V_cg)
    J = derivative(F, phi, dphi) 
    # Constant in front of time derivative
    C = Constant(e_v/(rho_w * g))


    ### Variational form 
    F = (theta*C*phi_dot - dot(grad(theta), q) + theta*(w - v - m))*dx
    F = theta*phi_dot*dx
    
    bdf = BDF(F, self.phi, self.phi_dot, self.model) 
    quit()
    
    self.F1 =F1
    self.J = J
    self.phi = phi
    self.model = model
    self.q = q
    self.shift = shift
    
    self.F3 = F3
    self.J3 = J3
    self.dt = dt
    #self.dt = dt
    #self.ts = ts
    
    
  def step_be_manual(self, dt):
    
    dt = 100.0
    #self.dt.assign(dt)
    #solve(self.F3 == 0, self.phi, self.model.d_bcs, J = self.J3, solver_parameters = self.model.newton_params)

    # Perturbation to be solved for
    dphi = Function(self.V_cg)
    # Old solution
    phi_old = Function(self.V_cg)
    
    #dt = 10.0
    self.shift.assign(1.0 / dt)
    eps = 1.0
    tol = 1e-5
    i = 0
    maxiter = 25
    
    while eps > tol and i < maxiter:
        i += 1
        # Compute F(x)
        F_x = assemble(-self.F3)
        self.model.d_bcs[0].apply(F_x)
        print "F: ", np.abs(F_x.array()).max()
        quit()
        print np.linalg.norm(F_x.array())
        # Compute J(x)
        J_x = assemble(self.J)
        self.model.d_bcs[0].apply(J_x)
        # Solve for perturbation
        solve(J_x, dphi.vector(), F_x)
        # Store old solution
        phi_old.assign(self.phi)
        # Update phi
        self.phi.vector()[:] += dphi.vector().array()
        # Update phi dot
        self.phi_dot.vector()[:] = (self.phi.vector().array() - phi_old.vector().array()) / dt
        
        #self.phi_dot.vector()[:] =
    
    
    
     # Step PDE for phi forward by dt. No constraints.
  def step(self, dt):
    self.step_be_manual(dt)
    
    """# Step the ODE forward
    self.ts.setTime(0.0)
    self.ts.setMaxTime(dt)
    
    self.ts.solve(self.phi_v)
    
    if self.MPI_rank == 0:
      print('steps %d (%d rejected)'
            % (self.ts.getStepNumber(), self.ts.getStepRejections()))
    
    # Apply changes to vector
    self.phi.vector().apply("insert")
    # Update phi
    self.model.update_phi()"""
    

    
    
  # Step PDE for phi forward by dt. Constrain using SNES solver. 
  def step_constrained(self):
    # Solve for potential
    (i, converged) = self.phi_solver.solve()
    # Update phi
    self.model.update_phi()  
    
