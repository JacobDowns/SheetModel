from dolfin import *
from constants import *
from pylab import *

"""
Starts from the steady state from sim_sliding_steady.py. This is a trough
simulation where everything is fixed except conductivity which is varied by 
an order of magnitude. 
"""
matplotlib.rcParams.update({'font.size': 16})

# Model input directory
in_dir = "out_cond_sensitivity1/checkpoint/"

# Load mesh and create function spaces
mesh = Mesh("inputs_sliding_law/mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

pfo = Function(V_cg)

### Run the simulation

# Time step
dt = 60.0 * 60.0 * 4.0

#ks = linspace(5e-4, 5e-2, 20)
ks = linspace(5e-4, 5e-2, 40)

t_10 = []
t_20 = []
t_50 = []
o_10 = []
o_20 = []
o_50 = []

for i in range(len(ks)):
  
  print i
  
  # Get the pfo for all k's
  File(in_dir + "pfo_" + str(i) + ".xml") >> pfo
  
  #dolfin.plot(model.pfo, interactive = True)
  
    # Record the pressure at different distances from the margin
  t_10.append(pfo([10e3, 10e3]))
  t_20.append(pfo([20e3, 10e3]))
  t_50.append(pfo([50e3, 10e3]))
  o_10.append(pfo([10e3, 2.5e3]))
  o_20.append(pfo([20e3, 2.5e3]))
  o_50.append(pfo([50e3, 2.5e3]))
  

title('Pressure as a Fraction of Overburden vs. Hydraulic Conductivity')
xlabel('Hydraulic Conductivity 'r'($m^{7/4}kg^{-1/2}$)')
ylabel('pfo')
xlim([0.0, max(ks)])
plot(ks, t_10, 'ro-', label = '10 km (Trough)')
plot(ks, t_20, 'bo-', label = '20 km (Trough)')
plot(ks, t_50, 'go-', label = '50 km (Trough)')
plot(ks, o_10, 'ro--', label = '10 km (Side)')
plot(ks, o_20, 'bo--', label = '20 km (Side)')
plot(ks, o_50, 'go--', label = '50 km (Side)')

legend()

show()