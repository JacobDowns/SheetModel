from dolfin import *
from dolfin import MPI, mpi_comm_world
from colored import fg, attr


class Solver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # Function space
    V_cg = model.V_cg
    # A reference to the model, which contains all the inputs we need
    self.model = model
    # Get melt rate
    m = model.m
    # Sheet height
    h = model.h
    # Previous sheet height
    h_prev = model.h_prev
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
    phi_reg = Constant(1e-15)

    

    ### Expressions used in the variational forms

    # Unknown 
    phi.assign(phi_0)
    # Expression for effective pressure in terms of potential
    N = phi_0 - phi
    # Flux vector
    q = -k * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(delta / 2.0) * grad(phi)
    # Opening term 
    #w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)
    w = u_b * (h_r - h) / Constant(l_r)    
    # Closing term
    v = Constant(A) * h * N**3
    # Time step
    dt = Constant(1.0)
    
    ### Variational forms
    
    # Test function
    theta = TestFunction(V_cg)
    
    # Variational form for phi
    F_phi = -dot(grad(theta), q) * dx + (w - v - m) * theta * dx
    d_phi = TrialFunction(V_cg)
    J_phi = derivative(F_phi, phi, d_phi) 

    # Setup the nonlinear problem for phi
    phi_problem = NonlinearVariationalProblem(F_phi, phi, model.d_bcs, J_phi)
    phi_problem.set_bounds(phi_min, phi_max)
    
    snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          "maximum_iterations": 250,
                                          "report": True,
                                          "error_on_nonconvergence": False, 
                                          "relative_tolerance" : 1e-6,
                                          "absolute_tolerance" : 1e-6}}
                                          
    phi_solver = NonlinearVariationalSolver(phi_problem)
    phi_solver.parameters.update(snes_solver_parameters)
    
    # Variational form for h
    F_h = ((h - h_prev) - (w - v) * dt) * theta * dx
    d_h = TrialFunction(V_cg)
    J_h = derivative(F_h, h, d_h) 
    
    
    ### Set local vars
    
    self.F_h = F_h
    self.J_h = J_h
    self.h = h
    self.dt = dt
    self.phi_solver = phi_solver
    self.q = q
    

  # Step PDE for phi forward by dt. Constrain using SNES solver. 
  def step_phi(self):
    # Solve for potential
    (i, converged) = self.phi_solver.solve()
    # Update phi
    self.model.update_phi()  
    
  
  # Step ODE for h forward by dt
  def step_h(self, dt):
    self.dt.assign(dt)
    # Solve ODE
    solve(self.F_h == 0, self.h, J = self.J_h, solver_parameters = self.model.newton_params)
    # Update h
    self.model.update_h()
    

  # Steps the potential forward with h fixed
  def step(self, dt):
    self.step_phi()
    self.step_h(dt)
