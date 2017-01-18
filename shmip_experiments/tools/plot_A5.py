# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 15:28:10 2017

@author: jake
"""

from pylab import *

# Load Mauro's tuning data
data = loadtxt('tuning_A5.txt', delimiter=",")
xs1 = data[:,0]
Ns1 = data[:,1]
qs1 = data[:,4]

# Load model output
xs2 = loadtxt('xs.txt')
Ns2 = loadtxt('A5_N_mean.txt')
qs2 = loadtxt('A5_q_mean.txt')

plot(xs1, Ns1, 'ko-', linewidth = 2)
plot(xs2, Ns2, 'ro-', linewidth = 2)
show()


