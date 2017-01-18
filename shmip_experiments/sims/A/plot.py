# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:23:56 2017

@author: jake
"""

from sqrt_steady_view import *
from pylab import *

# Load hdf5 file
view = SqrtSteadyView('results_A5_b/steady_A5.hdf5')

increments = 250
xs = np.linspace(1, 100e3, increments)

# Compute width averaged effective pressure
N_int = view.width_integrate_N()
Ns = [N_int([x, 20e3]) for x in xs]
Ns = np.array(Ns) / 20e3   

# Compute width averaged flux
q_int = view.width_integrate_q()
qs = [q_int([x, 20e3]) for x in xs]
qs = np.array(qs) / 20e3   


    
    
# Load Mauro's tuning data
data = loadtxt('tuning_A5.txt', delimiter=",")
xs_m = data[:,0]
Ns_m = data[:,1]
qs_m = data[:,4]


plot(xs, Ns, 'ko-', linewidth = 2)
plot(xs_m, Ns_m, 'ro-', linewidth = 2)
savefig('tune.png')

#tools.write_nc()
  
  

