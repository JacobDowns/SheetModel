# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:25:08 2015

@author: jake
"""

from dolfin import *
from pylab import *
from constants import *

matplotlib.rcParams.update({'font.size': 16})

h1_out = File('lag1/high_day_lag/h.pvd')
h2_out = File('lag1/low_day_lag/h.pvd')
h3_out = File('lag1/high_week_lag/h.pvd')
h4_out = File('lag1/low_week_lag/h.pvd')
h5_out = File('lag1/high_month_lag/h.pvd')
h6_out = File('lag1/low_month_lag/h.pvd')

pfo1_out = File('lag1/high_day_lag/pfo.pvd')
pfo2_out = File('lag1/low_day_lag/pfo.pvd')
pfo3_out = File('lag1/high_week_lag/pfo.pvd')
pfo4_out = File('lag1/low_week_lag/pfo.pvd')
pfo5_out = File('lag1/high_month_lag/pfo.pvd')
pfo6_out = File('lag1/low_month_lag/pfo.pvd')

input_file1 = HDF5File(mpi_comm_world(), 'lag1/high_day_lag.hdf5', 'r') 
input_file2 = HDF5File(mpi_comm_world(), 'lag1/low_day_lag.hdf5', 'r')
input_file3 = HDF5File(mpi_comm_world(), 'lag1/high_week_lag.hdf5', 'r')  
input_file4 = HDF5File(mpi_comm_world(), 'lag1/low_week_lag.hdf5', 'r')  
input_file5 = HDF5File(mpi_comm_world(), 'lag1/high_month_lag.hdf5', 'r')  
input_file6 = HDF5File(mpi_comm_world(), 'lag1/low_month_lag.hdf5', 'r')  

mesh = Mesh()
input_file1.read(mesh, "mesh", False) 
V_cg = FunctionSpace(mesh, "CG", 1)
H = Function(V_cg)
input_file1.read(H, 'H')
phi0 = project(pcs['rho_i'] * pcs['g'] * H, V_cg)

# Pressures for each smiulation
pfo1 = Function(V_cg)
pfo2 = Function(V_cg)
pfo3 = Function(V_cg)
pfo4 = Function(V_cg)
pfo5 = Function(V_cg)
pfo6 = Function(V_cg)

# Sheet heights for each simulation
h1 = Function(V_cg)
h2 = Function(V_cg)
h3 = Function(V_cg)
h4 = Function(V_cg)
h5 = Function(V_cg)
h6 = Function(V_cg)

# Conductivity for each simulation
k1 = Function(V_cg)
k2 = Function(V_cg)
k3 = Function(V_cg)
k4 = Function(V_cg)
k5 = Function(V_cg)
k6 = Function(V_cg)

# Points to record pressure
xs = [10e3, 20e3, 50e3]
ys = [10e3, 10e3, 10e3]

# List that contains lists of pressures for each point 
pfos1 = [[] for i in range(len(xs))]
pfos2 = [[] for i in range(len(xs))]
pfos3 = [[] for i in range(len(xs))]
pfos4 = [[] for i in range(len(xs))]
pfos5 = [[] for i in range(len(xs))]
pfos6 = [[] for i in range(len(xs))]

# List that contains lists of conductivities for each point
ks1 = [[] for i in range(len(xs))]
ks2 = [[] for i in range(len(xs))]
ks3 = [[] for i in range(len(xs))]
ks4 = [[] for i in range(len(xs))]
ks5 = [[] for i in range(len(xs))]
ks6 = [[] for i in range(len(xs))]

# Water volume over time
v1 = []
v2 = []
v3 = []
v4 = []
v5 = []
v6 = []

# Seconds per day
spd = pcs['spd']
# Seconds per year
spy = pcs['spy'] 
# Seconds per month
spm = pcs['spm']
# End time
T = 8.0 * spm
# Time step
dt = 60.0 * 60.0 * 8.0
# Times
ts = np.arange(0.0, T - spd, spd)

### Load all the data
for i in range(len(ts)):
  print ("i", i)
  
  # Load pressure
  input_file1.read(pfo1, "phi/vector_" + str(i))
  pfo1.vector()[:] = pfo1.vector().array() /phi0.vector().array()
  input_file2.read(pfo2, "phi/vector_" + str(i))
  pfo2.vector()[:] = pfo2.vector().array() /phi0.vector().array()
  input_file3.read(pfo3, "phi/vector_" + str(i))
  pfo3.vector()[:] = pfo3.vector().array() /phi0.vector().array()
  input_file4.read(pfo4, "phi/vector_" + str(i))
  pfo4.vector()[:] = pfo4.vector().array() /phi0.vector().array()
  input_file5.read(pfo5, "phi/vector_" + str(i))
  pfo5.vector()[:] = pfo5.vector().array() /phi0.vector().array()
  input_file6.read(pfo6, "phi/vector_" + str(i))
  pfo6.vector()[:] = pfo6.vector().array() /phi0.vector().array()
  
  # Load sheet height
  input_file1.read(h1, "h/vector_" + str(i))
  input_file2.read(h2, "h/vector_" + str(i))
  input_file3.read(h3, "h/vector_" + str(i))
  input_file4.read(h4, "h/vector_" + str(i))
  input_file5.read(h5, "h/vector_" + str(i))
  input_file6.read(h6, "h/vector_" + str(i))
  
  h1_out << h1
  h2_out << h2
  h3_out << h3
  h4_out << h4
  h5_out << h5
  h6_out << h6
   
    
  
  # Load sliding speed 
  #File(in_dir1 + "checkpoint/u_b_" + str(i) + ".xml") >> u_b
  
  for i in range(len(xs)):
    x = xs[i]
    y = ys[i]
    
    # Record pressures
    pfos1[i].append(pfo1([x, y]))
    pfos2[i].append(pfo2([x, y]))
    pfos3[i].append(pfo3([x, y]))
    pfos4[i].append(pfo4([x, y]))
    pfos5[i].append(pfo5([x, y]))
    pfos6[i].append(pfo6([x, y]))
  
  # Compute the total water volume 
  v1.append(assemble(h1 * dx))
  v2.append(assemble(h2 * dx))
  v3.append(assemble(h3 * dx))
  v4.append(assemble(h4 * dx))
  v5.append(assemble(h5 * dx))
  v6.append(assemble(h6 * dx))
 
# Plots all points for a given input array field
def plot_points(data, n = len(xs), l = None, s = None):
  for i in range(n):
    d = data[i]
    #plot(ts / spm, d, s[i], label = l[i], linewidth = 2)  
    plot(ts / spm, d, linewidth = 2)  
    
"""
savetxt("lag/pfos1.txt", array(pfos1))
savetxt("lag/pfos2.txt", array(pfos2))
savetxt("lag/pfos3.txt", array(pfos3))
savetxt("lag/pfos4.txt", array(pfos4))
savetxt("lag/pfos5.txt", array(pfos5))
savetxt("lag/pfos6.txt", array(pfos6))
  
savetxt("lag/k1.txt", array(ks1))
savetxt("lag/k2.txt", array(ks2))
savetxt("lag/k3.txt", array(ks3))
savetxt("lag/k4.txt", array(ks4))
savetxt("lag/k5.txt", array(ks5))
savetxt("lag/k6.txt", array(ks6))

savetxt("lag/v1.txt", array(v1))
savetxt("lag/v2.txt", array(v2))
savetxt("lag/v3.txt", array(v3))
savetxt("lag/v4.txt", array(v4))
savetxt("lag/v5.txt", array(v5))
savetxt("lag/v6.txt", array(v6))"""


### 1 day lag sheet volume
figure(0, (24, 7))
subplot(1,3,1)

title('Total Sheet Volume: Day Lag')
xlabel("Time (months)")
ylabel(r'Volume (m$^3)$')
xlim([0.0, max(ts / spm)])
#ylim([0.1e8, 1e8])
plot(ts / spm, v1, 'k', linewidth = 2, label = "High Melt")
plot(ts / spm, v2, 'k--', linewidth = 2, label = "Low Melt")
#plot_points(pfos1)
#plot_points(pfos3)
#plot_points(pfos5)
#legend(loc = 3)
grid(True)
#show()
#quit()


### 1 week lag sheet volume

subplot(1,3,2)

title('Total Sheet Volume: Week Lag')
xlabel("Time (months)")
#ylabel(r'Volume m$^3$')
xlim([0.0, max(ts / spm)])
#ylim([0.1e8, 1e8])
plot(ts / spm, v3, 'k', linewidth = 2)
plot(ts / spm, v4, 'k--', linewidth = 2)
grid(True)


### 1 month lag sheet volume

subplot(1,3,3)

title('Total Sheet Volume: Month Lag')
xlabel("Time (months)")
xlim([0.0, max(ts / spm)])
#ylim([0.1e8, 1e8])
plot(ts / spm, v5, 'k', linewidth = 2)
plot(ts / spm, v6, 'k--', linewidth = 2)
grid(True)

tight_layout()
show()
