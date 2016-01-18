# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:01:54 2015

@author: jake
"""

from dolfin import *
from constants import *

# Model input directory
in_dir = "inputs_sliding_law/"
# Output directory
out_dir = "out_sliding_law_steady/"
# Checkpoint directory
check_dir = out_dir + "checkpoint/"
# Process number
MPI_rank = MPI.rank(mpi_comm_world())

# Load mesh and create function spaces
mesh = Mesh(in_dir + "mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

phi = Function(V_cg)
File(check_dir + "phi_269.xml") >> phi

phi_0 = Function(V_cg)
File(in_dir + "phi_0.xml") >> phi_0

B = Function(V_cg)
File(in_dir + "B.xml") >> B

H = Function(V_cg)
File(in_dir + "H.xml") >> H

S = project(B + H, V_cg)

# Driving stress
rho_i = pcs['rho_i']
g = pcs['g'] 
tau_b = project(rho_i * g * H * sqrt(dot(grad(S), grad(S))))


N = project(phi_0 - phi, V_cg)

plot(N, interactive = True)


u_b = Function(V_cg)
u_b.vector()[:] = 2.5e8 * (tau_b.vector().array() / N.vector().array()**2)
#u_b_vals = u_b.vector().array()
#indexes_over = u_b.vector().array() > 100.0
#u_b.vector()[indexes_over] = 100.0

#u_b = project(A_s * (tau_b / N**2), V_cg)

plot(u_b, interactive = True)
#u_b = project(A_s
#plot(N , interactive = True)

"""
from pylab import *
from constants import * 

spy = pcs['spy']
def f(ts):
  return maximum(zeros(len(ts)), cos(((2*pi) / spy) * ts))
  
ts = linspace(0, spy, 250)
ys = f(ts)

plot(ts, ys)
show()"""