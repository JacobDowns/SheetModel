# -*- coding: utf-8 -*-
"""
Plot results.
"""

from pylab import *


# Load Mauro's tuning data
data = loadtxt('tuning_A3.txt', delimiter=",")
xs1 = data[:,0]
Ns1 = data[:,1]
qs1 = data[:,4]

# Load model output
xs2 = loadtxt('xs.txt')
Ns2 = loadtxt('N_mean.txt')
qs2 = loadtxt('q_mean.txt')

plot(xs1, qs1, 'ko-', linewidth = 2)
plot(xs2, qs2, 'ro-', linewidth = 2)
show()
