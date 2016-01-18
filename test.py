# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:57:49 2015

@author: jake
"""

"""
from dolfin import *
#from dolfin_adjoint import *

in_dir = "inputs_slope/"
mesh = Mesh(in_dir + "mesh.xml") 
V = FunctionSpace(mesh, "CG", 1)

u0 = Expression('sin(2.0)')

print type(u0)

if hasattr(u0, 't'):
 print "has t"

if isinstance(u0, dolfin.Expression):
  print "stuff"

plot(project(u0, V), interactive = True)

print type(u0)
quit()

for t in linspace(0.0, 1e6, 1):
  u0.t  = t
  
""" 

from pylab import *

spy = 60.0 * 60.0 * 24.0 * 365.0
ts = linspace(0.0, spy, 1000)
ms = 2.0 - cos( ((2.0 * pi) / spy) * ts)

plot(ts, ms, 'ro-')
show()