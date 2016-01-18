from dolfin import *
from constants import *
from pylab import *

matplotlib.rcParams.update({'font.size': 16})

# Model input directory
in_dir1 = "out_sliding_sensitivity1/checkpoint/"
in_dir2 = "out_cond_sensitivity1/checkpoint/"

# Load mesh and create function spaces
mesh = Mesh("inputs_sliding_law/mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

pfo1 = Function(V_cg)
pfo2 = Function(V_cg)

### Run the simulation

# Time step
dt = 60.0 * 60.0 * 4.0

#ks = linspace(5e-4, 5e-2, 20)
ubs = np.linspace(1.0, 100.0, 40)

t1_10 = []
t1_20 = []
t1_50 = []
o1_10 = []
o1_20 = []
o1_50 = []

t2_10 = []
t2_20 = []
t2_50 = []
o2_10 = []
o2_20 = []
o2_50 = []

for i in range(len(ubs)):
  
  print i
  
  # Get the pfo for all k's
  File(in_dir1 + "pfo_" + str(i) + ".xml") >> pfo1
  File(in_dir2 + "pfo_" + str(i) + ".xml") >> pfo2
  
  #dolfin.plot(model.pfo, interactive = True)
  
    # Record the pressure at different distances from the margin
  t1_10.append(pfo1([10e3, 10e3]))
  t1_20.append(pfo1([20e3, 10e3]))
  t1_50.append(pfo1([50e3, 10e3]))
  o1_10.append(pfo1([10e3, 2.5e3]))
  o1_20.append(pfo1([20e3, 2.5e3]))
  o1_50.append(pfo1([50e3, 2.5e3]))
  
  t2_10.append(pfo2([10e3, 10e3]))
  t2_20.append(pfo2([20e3, 10e3]))
  t2_50.append(pfo2([50e3, 10e3]))
  o2_10.append(pfo2([10e3, 2.5e3]))
  o2_20.append(pfo2([20e3, 2.5e3]))
  o2_50.append(pfo2([50e3, 2.5e3]))
  
ks = linspace(5e-4, 5e-2, 40)

print len(ks)
print 

subplot(1,2,1)
title('Pressure as a Fraction of Overburden vs. Hydraulic Conductivity')
xlabel('Hydraulic Conductivity 'r'($m^{7/4}kg^{-1/2}$)')
ylabel('pfo')
xlim([0.0, max(ks)])
plot(ks, t2_10, 'ro-', label = '10 km (Trough)')
plot(ks, t2_20, 'bo-', label = '20 km (Trough)')
plot(ks, t2_50, 'go-', label = '50 km (Trough)')
plot(ks, o2_10, 'ro--', label = '10 km (Side)')
plot(ks, o2_20, 'bo--', label = '20 km (Side)')
plot(ks, o2_50, 'go--', label = '50 km (Side)')
ylim([0.55, 1.05])

legend(loc = 3)


subplot(1,2,2)
title('Pressure as a Fraction of Overburden vs. Sliding Speed')
xlabel('Sliding Speed (m/a)$)')
ylabel('pfo')
xlim([0.0, max(ubs)])
plot(ubs, t1_10, 'ro-', label = '10 km (Trough)')
plot(ubs, t1_20, 'bo-', label = '20 km (Trough)')
plot(ubs, t1_50, 'go-', label = '50 km (Trough)')
plot(ubs, o1_10, 'ro--', label = '10 km (Side)')
plot(ubs, o1_20, 'bo--', label = '20 km (Side)')
plot(ubs, o1_50, 'go--', label = '50 km (Side)')
ylim([0.55, 1.05])

legend(loc = 3)

show()