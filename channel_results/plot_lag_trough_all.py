# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:25:08 2015

@author: jake
"""

from dolfin import *
from cr_tools import *
from constants import *
import pylab as plt


#matplotlib.rcParams.update({'font.size': 16})


plot_rect = PlotRect('lag_trough/high_no_lag.hdf5')
print plot_rect.get_total_m(0)
plot(plot_rect.get_phi(100), interactive = True)
Qs1 = plot_rect.get_int_channel_discharge(100)
qs1 = plot_rect.get_int_sheet_discharge(100)
plt.plot(plot_rect.xs, plt.array(Qs1) + plt.array(qs1))
plt.show()
quit()

input_file1 = HDF5File(mpi_comm_world(), 'lag_trough/high_day_lag.hdf5', 'r') 
input_file2 = HDF5File(mpi_comm_world(), 'lag_trough/low_day_lag.hdf5', 'r')
input_file3 = HDF5File(mpi_comm_world(), 'lag_trough/high_week_lag.hdf5', 'r')  
input_file4 = HDF5File(mpi_comm_world(), 'lag_trough/low_week_lag.hdf5', 'r')  
input_file5 = HDF5File(mpi_comm_world(), 'lag_trough/high_month_lag.hdf5', 'r')  
input_file6 = HDF5File(mpi_comm_world(), 'lag_trough/low_month_lag.hdf5', 'r')  

mesh = Mesh()
input_file1.read(mesh, "mesh", False) 
V_cg = FunctionSpace(mesh, "CG", 1)
V_cr = FunctionSpace(mesh, "CR", 1)
H = Function(V_cg)
B = Function(V_cg)
input_file1.read(H, 'H')
input_file1.read(B, 'B')


# Potential at 0 pressure
phi_m = project(pcs['rho_w'] * pcs['g'] * B, V_cg)
# Ice overburden pressure
p_i = project(pcs['rho_i'] * pcs['g'] * H, V_cg)
# Potential at overburden pressure
phi0 = project(phi_m + p_i, V_cg)


ff_S1 = FacetFunctionDouble(mesh)
ff_S2 = FacetFunctionDouble(mesh)
ff_S3 = FacetFunctionDouble(mesh)
ff_S4 = FacetFunctionDouble(mesh)
ff_S5 = FacetFunctionDouble(mesh)
ff_S6 = FacetFunctionDouble(mesh)

edge_lens = Function(V_cr)
input_file1.read(edge_lens, 'edge_lens')
cr_tools = CRTools(mesh, V_cg, V_cr, edge_lens)

# Pressures for each smiulation
pfo1 = Function(V_cg)
pfo2 = Function(V_cg)
pfo3 = Function(V_cg)
pfo4 = Function(V_cg)
pfo5 = Function(V_cg)
pfo6 = Function(V_cg)

# Channel height for each simulation
S1 = Function(V_cr)
S2 = Function(V_cr)
S3 = Function(V_cr)
S4 = Function(V_cr)
S5 = Function(V_cr)
S6 = Function(V_cr)


# Conductivity each simulation
k1 = Function(V_cg)
k2 = Function(V_cg)
k3 = Function(V_cg)
k4 = Function(V_cg)
k5 = Function(V_cg)
k6 = Function(V_cg)

# Sheet heights for each simulation
h1 = Function(V_cg)
h2 = Function(V_cg)
h3 = Function(V_cg)
h4 = Function(V_cg)
h5 = Function(V_cg)
h6 = Function(V_cg)

# Conductivity over time
ks1 = []
ks2 = []
ks3 = []
ks4 = []
ks5 = []
ks6 = []

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

# Channel volume
vS1 = []
vS2 = []
vS3 = []
vS4 = []
vS5 = []
vS6 = []

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
ts = np.arange(0.0, T - spd, dt)

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
  
  # Load channel height
  input_file1.read(S1, "S/vector_" + str(i))
  input_file2.read(S2, "S/vector_" + str(i))
  input_file3.read(S3, "S/vector_" + str(i))
  input_file4.read(S4, "S/vector_" + str(i))
  input_file5.read(S5, "S/vector_" + str(i))
  input_file6.read(S6, "S/vector_" + str(i))
  
  # Compute the amount of water stored in channel
  vS1.append(sum(S1.vector().array() * edge_lens.vector().array()))
  vS2.append(sum(S2.vector().array() * edge_lens.vector().array()))
  vS3.append(sum(S3.vector().array() * edge_lens.vector().array()))
  vS4.append(sum(S4.vector().array() * edge_lens.vector().array()))
  vS5.append(sum(S5.vector().array() * edge_lens.vector().array()))
  vS6.append(sum(S6.vector().array() * edge_lens.vector().array()))
  
  # Load conductivity  
  input_file1.read(k1, "m/vector_" + str(i))
  input_file2.read(k2, "m/vector_" + str(i))
  input_file3.read(k3, "m/vector_" + str(i))
  input_file4.read(k4, "m/vector_" + str(i))
  input_file5.read(k5, "m/vector_" + str(i))
  input_file6.read(k6, "m/vector_" + str(i))
  
  ks1.append(k1([60e3, 10e3]))
  ks2.append(k2([60e3, 10e3]))
  ks3.append(k3([60e3, 10e3]))
  ks4.append(k4([60e3, 10e3]))
  ks5.append(k5([60e3, 10e3]))
  ks6.append(k6([60e3, 10e3]))
  
  # Compute the total water volume 
  v1.append(assemble(h1 * dx))
  v2.append(assemble(h2 * dx))
  v3.append(assemble(h3 * dx))
  v4.append(assemble(h4 * dx))
  v5.append(assemble(h5 * dx))
  v6.append(assemble(h6 * dx))
  
  """"
  pfo1_out << pfo1
  pfo2_out << pfo2
  pfo3_out << pfo3
  pfo4_out << pfo4
  pfo5_out << pfo5
  pfo6_out << pfo6
  
  h1_out << h1
  h2_out << h2
  h3_out << h3
  h4_out << h4
  h5_out << h5
  h6_out << h6
  
  cr_tools.copy_cr_to_facet(S1, ff_S1)
  cr_tools.copy_cr_to_facet(S2, ff_S2)
  cr_tools.copy_cr_to_facet(S3, ff_S3)
  cr_tools.copy_cr_to_facet(S4, ff_S4)
  cr_tools.copy_cr_to_facet(S5, ff_S5)
  cr_tools.copy_cr_to_facet(S6, ff_S6)
  
  S1_out << ff_S1
  S2_out << ff_S2
  S3_out << ff_S3
  S4_out << ff_S4
  S5_out << ff_S5
  S6_out << ff_S6"""


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
plot(ts/ spm, array(vS1) + array(v1), 'k', linewidth = 2)

print vS1

"""plot(ts / spm, v2, 'r--', linewidth = 2, label = 'Low Melt, Day Lag')
plot(ts / spm, v3, 'g', linewidth = 2, label = 'High Melt, Week Lag')
plot(ts / spm, v4, 'g--', linewidth = 2, label = 'Low Melt, Week Lag')
plot(ts / spm, v5, 'b', linewidth = 2, label = 'High Melt, Month Lag')
plot(ts / spm, v6, 'b--', linewidth = 2, label = 'Low Melt, Month Lag')"""

grid(True)
legend(loc = 3)

show()
