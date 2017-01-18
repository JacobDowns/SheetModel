# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:23:56 2017

@author: jake
"""

from time_view import *
from pylab import *

# Load hdf5 file
view = TimeView('results_D1/out.hdf5')

### Sampe N at points

Ns1 = np.zeros(view.num_steps)
Ns2 = np.zeros(view.num_steps)
Ns3 = np.zeros(view.num_steps)
ts = np.zeros(view.num_steps)

for i in range(view.num_steps):
  print i
  N = view.get_N(i)
  t = view.get_t(i)

  Ns1[i] = N([10.0, 10.00])
  Ns2[i] = N([50.0, 10.00])
  Ns3[i] = N([90.0, 10.00])
  
  ts[i] = t

for i in range(view.num_steps / 365):
  ns1 = Ns2[i*365:(i+1)*365]
  plot(ns1)
  print len(ns1)  
  
#plot(ts, Ns1, 'k')
#plot(ts, Ns2, 'r')
#plot(ts, Ns3, 'b')
savefig('thing.png')
  

  
  

