from dolfin import *
from dolfin import MPI, mpi_comm_world
from petsc4py import PETSc

"""Solves ODE for the sheet height h with phi fixed."""

class HSolver():

  def __init__(self, model):
    
    # Process number
    self.MPI_rank = MPI.rank(mpi_comm_world())
    
    ### Get a few fields and parameters from the model
    
    # Effective pressure
    self.N = model.N   
    # Sliding speed
    self.u_b = model.u_b
    # Rate factor
    self.A = model.pcs['A']
    # Distance between bumps
    self.l_r = model.pcs['l_r']
    # Bump height
    h_r = model.h_r
    # Bump height vector
    h_r_n = h_r.vector().array()
    # Shet height function
    self.h = model.h
    

    ### Set up the sheet height ODE

    # h as petsc4py vec (metaphysically linked to h)
    self.h_v = as_backend_type(self.h.vector()).vec()
    comm = model.mesh.mpi_comm().tompi4py() 


    # Object that encapsulates ODE object for PETSc
    class ODE(object):
      
      def __init__(self):
        #  Part of  opening term that doesn't depend on h
        self.v_o_0 = None
        # Part of closure term that doesn't depend on h
        self.v_c_0 = None
      
      def rhs(self, ts, t, h, dhdt):
        h_n = h.getArray()
        # Sheet opening term
        v_o_n = self.v_o_0 * (h_r_n - h_n) 
        # Ensure that the opening term is non-negative
        v_o_n[v_o_n < 0.0] = 0.0
        # Sheet closure term
        v_c_n = self.v_c_0 * h_n
        # Set right hand side
        dhdt.setArray(v_o_n - v_c_n) 

    # Create PETSc time stepping solver  
    ode = ODE()
    ode_solver = PETSc.TS().create(comm=comm)
    ode_solver.setType(ode_solver.Type.RK)
    ode_solver.setRHSFunction(ode.rhs)
    ode_solver.setTime(0.0)
    ode_solver.setInitialTimeStep(0.0, 1.0)
    ode_solver.setTolerances(atol=1e-10, rtol=1e-15)
    ode_solver.setMaxSteps(100)
    ode_solver.setExactFinalTime(ode_solver.ExactFinalTimeOption.MATCHSTEP)


    ### Set local variables    
    
    self.ode = ode
    self.ode_solver = ode_solver
    self.model = model
    

  # Step the gap height h forward by dt
  def step(self, dt):
    # Precompute parts of the rhs that don't depend on h, but do depend on other
    # time dependent fields
    self.ode.v_o_0 = self.u_b.vector().array() / self.l_r
    self.ode.v_c_0 = self.A * self.N.vector().array()**3
    # Step the ODE forward
    self.ode_solver.setTime(0.0)
    self.ode_solver.setDuration(dt)
    self.ode_solver.solve(self.h_v)
    print('steps %d (%d rejected)'
          % (self.ode_solver.getStepNumber(), self.ode_solver.getStepRejections()))
    # Apply changes to vector
    self.model.h.vector().apply("insert")
  
 
