from dolfin import *
from dolfin_adjoint import *
set_log_level(ERROR)

import moola

n = 64
mesh = UnitSquareMesh(n, n)

cf = CellFunction("bool", mesh)
subdomain = CompiledSubDomain('std::abs(x[0]-0.5) < 0.25 && std::abs(x[1]-0.5) < 0.25')
subdomain.mark(cf, True)
mesh = refine(mesh, cf)

V = FunctionSpace(mesh, "CG", 1)
W = FunctionSpace(mesh, "DG", 0)

f = interpolate(Expression("x[0]+x[1]"), W, name = "Control")
u = Function(V, name='State')
v = TestFunction(V)

F = (inner(grad(u), grad(v)) - f*v)*dx
print f
print u
bc = DirichletBC(V, 0.0, "on_boundary")
solve(F == 0, u, bc)

x = SpatialCoordinate(mesh)
d = 1/(2*pi**2)*sin(pi*x[0])*sin(pi*x[1]) # the desired temperature profile

alpha = Constant(1e-6)
J = Functional(f*dx)
control = Control(f)
rf = ReducedFunctional(J, control)

minimize(rf, method = "L-BFGS-B", tol = 2e-08, options = {"disp": True})

"""
problem = rf.moola_problem()
f_moola = moola.DolfinPrimalVector(f)
solver = moola.NewtonCG(problem, f_moola, options={'gtol': 1e-9, 'maxiter': 20, 'display': 3, 'ncg_hesstol': 0} )
                                                   
sol = solver.solve()
f_opt = sol['control'].data

plot(f_opt, interactive=True, title="f_opt")

#Define the expressions of the analytical solution

f_analytic = Expression("sin(pi*x[0])*sin(pi*x[1])")
u_analytic = Expression("1/(2*pi*pi)*sin(pi*x[0])*sin(pi*x[1])")

f.assign(f_opt)
solve(F == 0, u, bc)
control_error = errornorm(f_analytic, f_opt)
state_error = errornorm(u_analytic, u)
print "h(min):           %e." % mesh.hmin()
print "Error in state:   %e." % state_error
print "Error in control: %e." % control_error"""

