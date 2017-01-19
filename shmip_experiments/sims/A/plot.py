# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:23:56 2017

@author: jake
"""

from sqrt_steady_view import *
from pylab import *

# Load hdf5 file
view1 = SqrtSteadyView('results_A_8e-35/steady_A5.hdf5')
view2 = SqrtSteadyView('results_A_9e-35/steady_A5.hdf5')
view3 = SqrtSteadyView('results_A_1e-25/steady_A5.hdf5')


# Load Mauro's tuning data
data = loadtxt('tuning_A5.txt', delimiter=",")
xs = data[:,0]
Ns_m = data[:,1]
qs_m = data[:,4]


increments = 250
#xs = np.linspace(1, 100e3, increments)

# Compute width averaged effective pressure
N_int1 = view1.width_integrate_N()
Ns1 = [N_int1([x, 20e3]) for x in xs]
Ns1 = np.array(Ns1) / 20e3   

N_int2 = view2.width_integrate_N()
Ns2 = [N_int2([x, 20e3]) for x in xs]
Ns2 = np.array(Ns2) / 20e3   

N_int3 = view3.width_integrate_N()
Ns3 = [N_int3([x, 20e3]) for x in xs]
Ns3 = np.array(Ns3) / 20e3   

print sum((Ns1 - Ns_m)**1)
print sum((Ns2 - Ns_m)**1)
print sum((Ns3 - Ns_m)**1)


print "Asf"
    

"""
print average(Ns1)
print average(Ns2)
print average(Ns3)
print average(Ns_m)"""

plot(xs, Ns1, 'ro-', linewidth = 2, ms = 2)
plot(xs, Ns2, 'go-', linewidth = 2, ms = 2)
plot(xs, Ns3, 'bo-', linewidth = 2, ms = 2)
plot(xs, Ns_m, 'ko-', linewidth = 2)
savefig('tune.png')

#tools.write_nc()
  
  

