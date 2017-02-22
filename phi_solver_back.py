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


    ### Set up time dependent pressure PDE

    F = (theta*C*phi_dot - dot(grad(theta), q) + theta*(w - v - m))*dx
    
    # BDF time stepper
    bdf = BDF(F, phi, phi_dot, model)

    self.model = model
    self.bdf = bdf
    self.phi_dot = phi_dot
    self.phi = phi
    self.q = q
    
    
  # Step PDE for phi forward by dt. No constraints.
  def step(self, dt):
    self.bdf.step_d1(dt)
    self.model.update_phi()
     
    
