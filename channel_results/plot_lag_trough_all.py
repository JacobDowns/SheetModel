# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:25:08 2015

@author: jake
"""

from dolfin import *
from pylab import *
from constants import *

matplotlib.rcParams.update({'font.size': 16})

h1_out = File('lag_trough/high_day_lag/h.pvd')
h2_out = File('lag_trough/low_day_lag/h.pvd')
h3_out = File('lag_trough/high_week_lag/h.pvd')
h4_out = File('lag_trough/low_week_lag/h.pvd')
h5_out = File('lag_trough/high_month_lag/h.pvd')
h6_out = File('lag_trough/low_month_lag/h.pvd')

phi1_out = File('lag_trough/high_day_lag/phi.pvd')
phi2_out = File('lag_trough/low_day_lag/phi.pvd')
phi3_out = File('lag_trough/high_week_lag/phi.pvd')
phi4_out = File('lag_trough/low_week_lag/phi.pvd')
phi5_out = File('lag_trough/high_month_lag/phi.pvd')
phi6_out = File('lag_trough/low_month_lag/phi.pvd')

input_file1 = HDF5File(mpi_comm_world(), 'lag_trough/high_day_lag.hdf5', 'r') 
input_file2 = HDF5File(mpi_comm_world(), 'lag_trough/low_day_lag.hdf5', 'r')
input_file3 = HDF5File(mpi_comm_world(), 'lag_trough/high_week_lag.hdf5', 'r')  
input_file4 = HDF5File(mpi_comm_world(), 'lag_trough/low_week_lag.hdf5', 'r')  
input_file5 = HDF5File(mpi_comm_world(), 'lag_trough/high_month_lag.hdf5', 'r')  
input_file6 = HDF5File(mpi_comm_world(), 'lag_trough/low_month_lag.hdf5', 'r')  

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

# Average pressure over time
avg_pfos1 = []
avg_pfos2 = []
avg_pfos3 = []
avg_pfos4 = []
avg_pfos5 = []
avg_pfos6 = []
f = interpolate(Constant(1.0), V_cg)
area = assemble(f * dx)

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
  phi1_out << pfo1
  pfo1.vector()[:] = pfo1.vector().array() /phi0.vector().array()
  input_file2.read(pfo2, "phi/vector_" + str(i))
  phi2_out << pfo2
  pfo2.vector()[:] = pfo2.vector().array() /phi0.vector().array()
  input_file3.read(pfo3, "phi/vector_" + str(i))
  phi3_out << pfo3
  pfo3.vector()[:] = pfo3.vector().array() /phi0.vector().array()
  input_file4.read(pfo4, "phi/vector_" + str(i))
  phi4_out << pfo4
  pfo4.vector()[:] = pfo4.vector().array() /phi0.vector().array()
  input_file5.read(pfo5, "phi/vector_" + str(i))
  phi5_out << pfo5
  pfo5.vector()[:] = pfo5.vector().array() /phi0.vector().array()
  input_file6.read(pfo6, "phi/vector_" + str(i))
  phi6_out << pfo6
  pfo6.vector()[:] = pfo6.vector().array() /phi0.vector().array()
  
  avg_pfos1.append(assemble(pfo1 * dx) / area)
  avg_pfos2.append(assemble(pfo2 * dx) / area)
  avg_pfos3.append(assemble(pfo3 * dx) / area)
  avg_pfos4.append(assemble(pfo4 * dx) / area)
  avg_pfos5.append(assemble(pfo5 * dx) / area)
  avg_pfos6.append(assemble(pfo6 * dx) / area)
  
  # Load sheet height
  input_file1.read(h1, "h/vector_" + str(i))
  input_file2.read(h2, "h/vector_" + str(i))
  input_file3.read(h3, "h/vector_" + str(i))
  input_file4.read(h4, "h/vector_" + str(i))
  input_file5.read(h5, "h/vector_" + str(i))
  input_file6.read(h6, "h/vector_" + str(i))
  
  # Compute the total water volume 
  v1.append(assemble(h1 * dx))
  v2.append(assemble(h2 * dx))
  v3.append(assemble(h3 * dx))
  v4.append(assemble(h4 * dx))
  v5.append(assemble(h5 * dx))
  v6.append(assemble(h6 * dx))


figure(0, (24, 7))

### Plot average pressure over time
subplot(1,2,1)

title('Average PFO')
xlabel("Time (months)")
xlim([0.0, max(ts / spm)])

plot(ts / spm, avg_pfos1, 'r', linewidth = 2, label = 'High Melt, Day Lag')
plot(ts / spm, avg_pfos2, 'r--', linewidth = 2, label = 'Low Melt, Day Lag')
plot(ts / spm, avg_pfos3, 'g', linewidth = 2, label = 'High Melt, Week Lag')
plot(ts / spm, avg_pfos4, 'g--', linewidth = 2, label = 'Low Melt, Week Lag')
plot(ts / spm, avg_pfos5, 'b', linewidth = 2, label = 'High Melt, Month Lag')
plot(ts / spm, avg_pfos6, 'b--', linewidth = 2, label = 'Low Melt, Month Lag')

grid(True)
legend(loc = 3)


### Plot sheet volume over time

subplot(1,2,2)

title('Sheet Volume')
xlabel("Time (months)")
ylabel(r'Volume (m$^3$)')
xlim([0.0, max(ts / spm)])
ylim([2e7, 9e7])

plot(ts / spm, v1, 'r', linewidth = 2, label = 'High Melt, Day Lag')
plot(ts / spm, v2, 'r--', linewidth = 2, label = 'Low Melt, Day Lag')
plot(ts / spm, v3, 'g', linewidth = 2, label = 'High Melt, Week Lag')
plot(ts / spm, v4, 'g--', linewidth = 2, label = 'Low Melt, Week Lag')
plot(ts / spm, v5, 'b', linewidth = 2, label = 'High Melt, Month Lag')
plot(ts / spm, v6, 'b--', linewidth = 2, label = 'Low Melt, Month Lag')

grid(True)
legend(loc = 3)

show()
