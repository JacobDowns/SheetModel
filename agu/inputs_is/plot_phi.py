from dolfin import *

mesh = Mesh("mesh.xml")
V_cg = FunctionSpace(mesh, 'CG', 1)
V2 = VectorFunctionSpace(mesh, 'CG', 1)

pfo = Function(V_cg)
File("pfo_6.xml") >> pfo

P_i = Function(V_cg)
File("p_i.xml") >> P_i

phi_m = Function(V_cg)
File("phi_m.xml") >> phi_m

P_w = pfo * P_i

phi = phi_m + P_w

plot(phi, interactive = True)
plot(grad(phi), interactive = True)
File("flux.pvd") << project(-grad(phi), V2)
