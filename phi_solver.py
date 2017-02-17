from dolfin import *
from dolfin import MPI, mpi_comm_world
from petsc4py import PETSc

""" Solves phi with h fixed."""

class PhiSolver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # A reference to the model, which contains all the inputs we need
    self.model = model
    # Function space
    V_cg = model.V_cg
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
    self.phi1 = phi1
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

    # Unknown 
    #phi.assign(phi_0)
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


    ### Variational form used to compute f for PETSc
    F1 = (theta*C*phi_dot - dot(grad(theta), q) + theta*(w - v - m))*dx
    # f as fenics vector
    f_f = assemble(F1)
    # f as petsc vector
    f_p = as_backend_type(f_f).vec()
    
    
    ### Variational form used to compute Jacobian for PETSc
    shift = Constant(1.0)
    F2 = (theta*shift*C*phi - dot(grad(theta), q) + theta*(w - v - m))*dx

    # Jacobian as form expression
    dphi = TrialFunction(V_cg)
    J = derivative(F2, phi, dphi)
    # Jacobian as fenics vector
    J_f = assemble(J)
    # Jacobian as PETSc matrix
    J_p = as_backend_type(J_f).mat()
    
    
    ### Variational form for backward Euler
    dt = Constant(1.0)
    F3 = (theta*C*(phi - phi1) + dt*(-dot(grad(theta), q) + theta*(w - v - m)))*dx
    J_F3 = derivative(F3, phi, dphi) 
    
    self.dt = dt
    self.J_F3 = J_F3
    self.F3 = F3
    
    ### Set up 
    
    for bc in model.d_bcs:
      bc.apply(phi.vector())

    # Object that encapsulates ODE object for PETSc
    class ODE(object):

      def evalFunction(self, ts, t, x, xdot, f):
        #if MPI_rank == 0:
        #print t

        phi_dot.vector().set_local(xdot.getArray())
        phi_dot.vector().apply("insert")     
        
        phi.vector().set_local(x.getArray())
        phi.vector().apply("insert")
        
        assemble(F1, tensor = f_f)
        
        for bc in model.d_bcs:
          bc.apply(f_f)
          
        #print f.getArray().max()
       # print (f.getArray() - f_f.array()).max()
        #print 
      
      def evalJacobian(self, ts, t, x, xdot, a, A, B):

        shift.assign(a)
        
        phi_dot.vector().set_local(xdot.getArray())
        phi_dot.vector().apply("insert")  
        
        phi.vector().set_local(x.getArray())
        phi.vector().apply("insert")

        #for bc in model.d_bcs:
        #  bc.apply(phi.vector())
          
        assemble(J, tensor = J_f)
        
        for bc in model.d_bcs:
          bc.apply(J_f)
        
        # If operator is different from preconditioning matrix
        #print A.getValues([0,1], [0,1])
        if A != B: 
          A.assemble()
          # Same nonzero pattern
          return True 
    
    
    self.phi_v = as_backend_type(phi.vector()).vec()
    comm = model.mesh.mpi_comm().tompi4py() 
    
    ode = ODE()
    ts = PETSc.TS().create(comm=comm)
    ts.setProblemType(ts.ProblemType.NONLINEAR)
    ts.setEquationType(ts.EquationType.IMPLICIT)
    ts.setType(ts.Type.BDF)
    ts.setIFunction(ode.evalFunction, f_p)
    ts.setIJacobian(ode.evalJacobian, J_p)
    ts.setTime(0.0)
    ts.setInitialTimeStep(0.0, 60.0 * 10.0)
    #ts.setTimeStep(60.0)    
    ts.setTolerances(atol=1e-5, rtol=1e-9)
    ts.setMaxSteps(100)
    ts.setExactFinalTime(ts.ExactFinalTimeOption.MATCHSTEP)
    ts.setMaxSNESFailures(-1)
    
    snes = ts.getSNES()             # Nonlinear solver
    snes.setTolerances(max_it=10)   # Stop nonlinear solve after 10 iterations (TS will retry with shorter step)
    ksp = snes.getKSP()             # Linear solver
    ksp.setType(ksp.Type.CG)        # Conjugate gradients
    pc = ksp.getPC()                # Preconditioner
    if False:                       # Configure algebraic multigrid, could use run-time options instead
        pc.setType(pc.Type.GAMG)    # PETSc's native AMG implementation, mostly based on smoothed aggregation
        OptDB['mg_coarse_pc_type'] = 'svd' # more specific multigrid options
        OptDB['mg_levels_pc_type'] = 'sor'
    
    self.F = F
    self.J = J
    self.phi = phi
    self.model = model
    self.q = q
    self.ts = ts
    
    
  # Step PDE for phi forward by dt using backward Euler
  def step_be(self, dt):
    # Assign time step
    self.dt.assign(dt)
    # Solve for potential
    solve(self.F3 == 0, self.phi, self.model.d_bcs, J = self.J_F3, solver_parameters = self.model.newton_params)
    # Update phi1
    self.phi1.assign(self.phi)
    # Update fields derived from phi
    self.model.update_phi()
    
    
  # Step PDE for phi forward by dt. No constraints.
  def step(self, dt):
    self.step_be(dt)

    """"
    # Step the ODE forward
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
    
