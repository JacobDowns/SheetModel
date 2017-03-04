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
    # 6 x 6 identity
    self.I = identity(6)
    # Special matrix for computing n-th derivatives of Newton polynomial
    self.M = array([[ 0.,  0.,  0.,  0.,  0.,  0.],
                    [ 1.,  0.,  0.,  0.,  0.,  0.],
                    [ 0.,  2.,  0.,  0.,  0.,  0.],
                    [ 0.,  0.,  3.,  0.,  0.,  0.],
                    [ 0.,  0.,  0.,  4.,  0.,  0.],
                    [ 0.,  0.,  0.,  0.,  5.,  0.]])
      
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
    
    # Derivatives of N0 [dN0^(0), ..., dN0^(5)]
    dNi_k = zeros(k + 1)
    dNi_k[0] = 1.0
    
    # ci's
    coefs = self.get_dd_coefs_n(0, times)*dNi_k[k]

    # Compute derivatives of N_(i+1)
    for i in range(len(times)-1):
      dNi_k = dot(((t - times[i])*self.I[0:(k+1), 0:(k+1)] + self.M[0:(k+1), 0:(k+1)]), dNi_k)
      coefs += self.get_dd_coefs_n(i+1, times)*dNi_k[k]
    
    return coefs
    
  
    
      
from pylab import *
bdf_helper = NewtonPoly()


ts = linspace(2.5, 6.0, 100)
times = array([6.0, 5.0, 4.0, 3.0, 0.5])
ys = array([5.0, 1.0, 7.0, -1., 5.0])



p0 = []
p1 = []
p2 = []
for t in ts:
  p0.append(dot(bdf_helper.get_dpk_dt_coefs_at_t(0, times, t), ys))
  p1.append(dot(bdf_helper.get_dpk_dt_coefs_at_t(1, times, t), ys))
  p2.append(dot(bdf_helper.get_dpk_dt_coefs_at_t(2, times, t), ys))
  
plot(ts, 0.0*ts, 'k')
plot(ts, p0)
plot(ts, p1)
plot(ts, p2)
plot(times, ys, 'ko-')
xlim([ts.min(), ts.max()])
show()