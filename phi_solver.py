from dolfin import *
from dolfin import MPI, mpi_comm_world

""" Solves phi with h fixed."""

class PhiSolver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # A reference to the model, which contains all the inputs we need
    self.model = model
    # Get melt rate
    m = model.m
    # Sheet height
    h = model.h
    # Basal sliding speed
    u_b = model.u_b
    # Potential
    phi = model.phi
    self.phi = phi
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
    theta = TestFunction(model.V_cg)
    # Variational form for the PDE
    F = -dot(grad(theta), q) * dx + (w - v - m) * theta * dx
    # Get the Jacobian
    dphi = TrialFunction(model.V_cg)
    J = derivative(F, phi, dphi) 

    
    # Setup the nonlinear problem for the hydraulic potential
    phi_problem = NonlinearVariationalProblem(F, phi, model.d_bcs, J)
    phi_problem.set_bounds(phi_min, phi_max)
    
    snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          "maximum_iterations": 100,
                                          'line_search': 'basic',
                                          "report": True,
                                          "error_on_nonconvergence": False, 
                                          "relative_tolerance" : 1e-11,
                                          "absolute_tolerance" : 1e-9}}
                      
    # Set object variables                  
    self.phi_solver = NonlinearVariationalSolver(phi_problem)
    self.phi_solver.parameters.update(snes_solver_parameters)
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
    
