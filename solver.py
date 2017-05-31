from dolfin import *
from dolfin import MPI, mpi_comm_world
import numpy as np
import ufl
from petsc4py import PETSc

""" Steps forward (phi, h)"""
ufl.algorithms.apply_derivatives.CONDITIONAL_WORKAROUND = True

class Solver(object):
  
  def __init__(self, model):
    
    # Process number    
    self.MPI_rank = MPI.rank(mpi_comm_world())
    # A reference to the model, which contains all the inputs we need
    self.model = model
    # MPI comm object    
    comm = model.mesh.mpi_comm().tompi4py() 
    
    
    ### Get all the variables, inputs, and constants we need from the model
    
    # Mixed function space
    V = model.V
    # Combined (phi, h) unknown
    U = model.U
    # CG Function space
    V_cg = self.model.V_cg
    # Hydraulic potential 
    phi, h = split(U)
    # Derivative of combined unknown
    U_dot = Function(self.model.V)
    # Time derivative of h
    h_dot = Function(V_cg)
    # Time derivative of phi
    phi_dot = Function(V_cg)
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
    # Apply boundary conditions for pressure
    [bc.apply(U.vector()) for bc in model.d_bcs]
    # Vector of error tolerances    

    # Expression for effective pressure in terms of potential
    N = phi_0 - phi
    # Flux vector
    q = -k * h**alpha * (dot(grad(phi), grad(phi)) + phi_reg)**(delta / 2.0) * grad(phi)
    # Opening term 
    w = conditional(gt(h_r - h, 0.0), u_b * (h_r - h) / Constant(l_r), 0.0)
    # Closing term
    v = Constant(A) * h * N**3
    
    
    ### Form for DAE
    
    # Test functions
    theta1, theta2 = TestFunctions(model.V)
    # PDE part of F
    F_form = -dot(grad(theta1), q)*dx + (w - v - m)*theta1*dx
    # ODE part of F
    F_form += (h_dot - w + v)*theta2*dx
    
        
    ### Form for Jacobian
        
    # Get the Jacobian
    dU = TrialFunction(model.V)
    # Note: the following trick of substituting h_dot by h works because the 
    # Jacobian wrt to h_dot is constant
    shift = Constant(1.0)
    J_form = -dot(grad(theta1), q)*dx + (w - v - m)*theta1*dx
    J_form += (shift*h - w + v)*theta2*dx
    J = derivative(J_form, U, dU)
    

    ### Setup time stepper 

    # Create PETSc vectors of atols and rtols
    dif_atol = 1e-6
    alg_atol = 1e16
    dif_rtol = 1e-8
    alg_rtol = 1e16
    
    dif_atol_cg = project(Constant(dif_atol), V_cg)
    dif_rtol_cg = project(Constant(dif_rtol), V_cg)
    alg_atol_cg = project(Constant(alg_atol), V_cg)
    alg_rtol_cg = project(Constant(alg_rtol), V_cg)  
    
    atol = Function(V)
    rtol = Function(V)
    model.assign_to_mixed.assign(atol, [alg_atol_cg, dif_atol_cg])
    model.assign_to_mixed.assign(rtol, [alg_rtol_cg, dif_rtol_cg])
    atol_p = as_backend_type(atol.vector()).vec()
    rtol_p = as_backend_type(rtol.vector()).vec()
    
    # F as a fenics vector
    F_f = assemble(F_form)
    [bc.apply(F_f) for bc in model.d_bcs]
    F_p = as_backend_type(F_f).vec()
    F_p.assemble()
    J_f = assemble(J)
    [bc.apply(J_f) for bc in model.d_bcs]
    #print U.vector().array()
    #print J_f.array()
    #print
    J_p = as_backend_type(J_f).mat()
    J_p.assemble()
    
    # U as a petsc vector
    U_p = as_backend_type(U.vector()).vec().duplicate()
    U_p.setArray(U.vector().array())
    U_p.assemble()
    
    # ODE object for PETSc
    class ODE(object):

      def evalFunction(self, ts, t, u, u_dot, f):   
        print t
        # Assign h_dot
        U_dot.vector().set_local(u_dot.getArray())
        U_dot.vector().apply("insert")        
        model.assign_from_mixed.assign([phi_dot, h_dot], U_dot)
        
        # Assign U
        U.vector().set_local(u.getArray())
        U.vector().apply("insert")
        
        # Assemble F
        assemble(F_form, tensor = F_f)
        [bc.apply(F_f) for bc in model.d_bcs]
        f.setArray(as_backend_type(F_f).vec().getArray())
        print ("fmax", f.getArray().max())
        
      
      def evalJacobian(self, ts, t, u, u_dot, a, A, B): 
        # Assign U, no need to assign U_dot because derivatives wrt U_dot are constant       
        U.vector().set_local(u.getArray())
        U.vector().apply("insert")
        # Assign shift
        shift.assign(a)
        
        # Assemble J
        assemble(J, tensor = J_f)
        [bc.apply(J_f) for bc in model.d_bcs]
        
        print np.max(as_backend_type(J_f).mat().getValuesCSR()[2] - B.getValuesCSR()[2])
        print
        
        if A != B: 
          A.assemble()
      
          
    ode = ODE()
    ts = PETSc.TS().create(comm=comm)
    ts.setProblemType(ts.ProblemType.NONLINEAR)
    ts.setEquationType(ts.EquationType.IMPLICIT)
    ts.setType(ts.Type.ROSW)
    ts.setIFunction(ode.evalFunction, F_p)
    ts.setIJacobian(ode.evalJacobian, J_p)
    ts.setTime(0.0)
    ts.setInitialTimeStep(0.0, 60.0)   
    ts.setTolerances(vatol = atol_p, vrtol = rtol_p)
    ts.setMaxSteps(500)
    ts.setExactFinalTime(ts.ExactFinalTimeOption.MATCHSTEP)
    ts.setMaxSNESFailures(-1)
    
    snes = ts.getSNES()             # Nonlinear solver
    snes.setTolerances(rtol = 1e-4, atol=1e-4, max_it=20)   # Stop nonlinear solve after 10 iterations (TS will retry with shorter step)
    ksp = snes.getKSP()             # Linear solver
    
                      
    # Set object variables                  
    self.U = U
    self.model = model
    self.q = q
    self.ts = ts
    self.U_p = U_p
    self.U_dot = U_dot
    
    
   # Step PDE for phi forward by dt. No constraints.
  def step(self, dt):
    # Step the ODE forward
    self.ts.setTime(0.0)
    self.ts.setMaxTime(dt)
    self.ts.solve(self.U_p)
    
    if self.MPI_rank == 0:
      print('steps %d (%d rejected)'
            % (self.ts.getStepNumber(), self.ts.getStepRejections()))
    
    # Apply changes to vector
    self.U.vector().set_local(self.U_p.getArray())
    self.U.vector().apply("insert")
    # Update phi
    self.model.update_phi()
