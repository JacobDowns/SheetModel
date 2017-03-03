# -*- coding: utf-8 -*-
"""
This class has a bunch of functions for dealing with Newton interpolating 
polynomials. A Newton polynomial passes through the points (t0, y0) ... (tk, yk).
A Newton polynomial can be written in terms of divided differences as:

  p(t) = [y0] + [y0, y1] * N_1(t) + ... + [y0, ..., yk] * N_k(t)
  
where N_1(t) = (t - t0), N_2(t) = (t - t0)(t - t1) and so on. 
"""

from numpy import *
from scipy.ndimage.interpolation import shift

class NewtonPoly(object):
  
  def __init__(self):
    pass
  
  """
  If p is the Newton polynomial passing through y0, ... yk at the given times
  t0 ... tk then we can evaluate p(t) as a linear combination 
  c0*y0 + ... + ck*yk. This function returns the coefficients ci. 
  """
  def get_p_coefs_at_t(self, times, t):
    coefs = self.get_dd_coefs_n(len(times)-1, times)

    for j in range(len(times)-1, 0,-1):
      coefs = self.get_dd_coefs_n(j-1, times) + (t - times[j-1])*coefs
      
    return coefs
    
  """
  If p is the Newton polynomial passing through y0, ... yk at the given times
  t0 ... tk then we can evaluate p'(t) as a linear combination 
  c0*y0 + ... + ck*yk. This function returns the coefficients ci. 
  """
  def get_dp_dt_coefs_at_t(self, times, t):
    coefs_p = self.get_dd_coefs_n(len(times)-1, times)
    coefs_pp = zeros(len(times))
    
    for j in range(len(times)-1, 0,-1):
      coefs_pp = coefs_p + (t - times[j-1])*coefs_pp      
      coefs_p = self.get_dd_coefs_n(j-1, times) + (t - times[j-1])*coefs_p
      
    return coefs_pp
      
  """ 
  A divided difference [y0, ... , yn] returns a linear comibination
  of c0*y0, ..., cn*yn. This function returns the coefficients ci.
  """
  def get_dd_coefs_n(self, n, times):
    n = min(n, len(times))
    weights = zeros(len(times))
    
    for i in range(n+1):
      wi = (times[i] - times[:n+1])
      wi[i] = 1.0
      wi = prod(wi)
      weights[i] = 1.0 / wi
      
    return weights
    
    
  """
  If p is the Newton polynomial passing through y0, ... yk at the given times
  t0 ... tk then we can evaluate p^(n)(t) as a linear combination 
  c0*y0 + ... + ck*yk. This function returns the coefficients ci. 
  """
  
  def get_dpk_dt_coefs_at_t(self, k, times, t):
    
    k = min(k, len(times))
    
    # First we'll build an array that stores [dN_k^(k), ..., dN_n^(k)]
    dN_dts = zeros(len(times) - k)
     
    # ith derivative of Ni 
    dNi_i1_dt = 1.0
    # (i-1)th derivative of Ni
    dNi_i0_dt = 0.0
    
    
    for i in range(1,k+1):
      dNi_i0_dt = (t -times[i - 1])*dNi_i1_dt + (i - 1.0)*dNi_i0_dt
      dNi_i1_dt = i*dNi_i1_dt

    dN_dts[0] = dNi_i1_dt
    
    print len(times) - k
    weights[0] = dNi_i1_dt
    quit()
    
    #for i in range(k, len(times)):
      
    
      # N_(k-1)^(k-1) is first non-zero derivative 
    
  
    
      

bdf_helper = NewtonPoly()
times = array([6.0, 5.0, 4.0, 3.0])

bdf_helper.get_dpk_dt_coefs_at_t(1, times, 6.1)
quit()
from pylab import *

ys = array([5.0, 1.0, 7.0, 0.5])


ts = linspace(2.5, 6.0, 100)
ps = []
pps = []
for t in ts:
  ps.append(dot(bdf_helper.get_p_coefs_at_t(times, t), ys))
  pps.append(dot(bdf_helper.get_dp_dt_coefs_at_t(times, t), ys))
  
plot(ts, ps)
plot(ts, pps)
plot(times, ys, 'ko-')
xlim([ts.min(), ts.max()])
show()
