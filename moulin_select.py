# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 17:07:46 2015

@author: jake
"""

from dolfin import *
from numpy import *

in_dir = "inputs_trough1/"
mesh = Mesh(in_dir + "mesh_60_20.xml")
V_cg = FunctionSpace(mesh, "CG", 1)


### First select som points where we'd like moulins

xs = arange(5000, 60000, 5000.0)
ys = random.rand(len(xs)) * 20000.0


### Next find create a list of dofs that are close to the desired moulins

vertex_coords = V_cg.dofmap().tabulate_all_coordinates(mesh).reshape(V_cg.dim(), 2)
dofs = []

for i in range(len(xs)):
  x = xs[i]
  y = ys[i]
  
  dist = vertex_coords - array([x,y])
  dist = dist**2
  dist = sum(dist, axis = 1)
  
  # Get the index of the vertex closest to the desired point
  vertex_index = dist.argmin()
  dofs.append(vertex_index)

print dofs

### Calculate the desired melt rate at each moulin 

# Melt function
def m(x):
  return (1.0 + (4.0 * (60000.0 - x) / 60000.0)) / 31536000.0

# Array with the starting and stopping x coordinates of each moulin domain

domain_starts = array(xs)
domain_starts[0] = 0.0
domain_starts = domain_starts[:]

domain_stops = array(xs)
domain_stops = domain_stops[1:]
domain_stops = append(domain_stops, 60000.0)

print domain_starts
print domain_stops


domain_lens = domain_stops - domain_starts

melt_starts = m(domain_starts)
melt_stops = m(domain_stops)

w = 20000.0
moulin_melts = []
for i in range(len(domain_starts)):
  l = domain_lens[i]  
  m1 = melt_starts[i]
  m2 = melt_stops[i]
  
  total_melt = l * w * m2
  total_melt += l * w * ((m1 - m2) / 2.0)
  
  moulin_melts.append(total_melt)

print sum(moulin_melts) 

### Compute base integral under each hat function
m = Function(V_cg)

hat_integrals = []
for dof in dofs:
  m.vector()[:] = zeros(V_cg.dim())  
  m.vector()[dof] = 1.0
  hat_integrals.append(assemble(m * dx))

print len(hat_integrals)
# From the base integrals we can computes the desired dof values
dof_vals = array(moulin_melts) / array(hat_integrals)

print dof_vals

m.vector()[dofs] = dof_vals
print assemble(m * dx)

plot(m, interactive = True)

File("inputs_trough1/m_point2.pvd") << m
File("inputs_trough1/m_point2.xml") << m
  