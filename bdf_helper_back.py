# -*- coding: utf-8 -*-
"""
Class to compute time dependent coefficients and truncation error for a k-th order 
BDF method (k <= 5).

phi(t) is the interpolating polynomial that passes through 
(t0, y0), (t1, y1), ... , (t5, y5).
phidot(t) is the derivative of the interpolating polynomial wrt t. It can be written
in terms of divided differences as
phidot(t) = [y0, y1] + [y0, y1, y2](t0-t1) + [y0, y1, y2, y3](t0-t1)(t0-t2) + ...
which can be rewritten as a linear cominbation of y0 ... yn
"""

from numpy import *
from scipy.ndimage.interpolation import shift

class BDF_Helper(object):
  
  def __init__(self):
    # Times for last 6 solutions   
    #[t0, t1, t2, t3, t4, t5]
    self.ts = zeros(6)
    
  
  # Suppose p(t) is the interpolating polynomial passing through 
  # (t0, y0) ... (t_(k-1), y_(k-1) p evaluated at t can be expressed as a 
  # linear combination a0*y0 + ... + a_(k-1)*y_(k-1).This function returns 
  # the coefficients ai
  def get_pred_coefs_at_t(self, k, times, t):
    coefs = zeros(6)
    
    for i in range(k):
      n_i = (t - times[0:k])
      n_i[i] = 1.0

      d_i = (times[i] - times[0:k])
      d_i[i] = 1.0
        
      coefs[i] = prod(n_i / d_i)
    
    return coefs
    
  """ 
  -k is degree
  -prev_times = [t0, t1, ... , t5] is an array of times of previous solutions
  
  Let p(t) be the so called "predictor" polynomial that passes through the points
  (t0, y0), ... (t_(k-1), y_(k-1))
  
  The "corrector" polynomial c(t) is such that c(t0) = y0, and c(t0 - kh) = p(t0 - kh)
  for k = [0, ..., k-1]. That is it agrees with p(t) on a regular grid. For a BDF
  method, we want to enforce (dc/dt)(t0+h) = f(y0, t0+h). 
  
  (dc/dt)(t0+h) can be expressed as a linear combination 
  
    (dc/dt)(t0+h) = a0*y0 + a1*p(t0) + a2*p(t0-h) + a3*p(t0-2h) + ... 
  
  In turn p(t0), p(t0-h) can be written as linear combinations os y1, y2, ...
  so that we have 
  
    (dc/dt)(t0+h) = a0*y0 + a1'*y1 + a2'*y2 + ...
  
  This function returns the coefficients a0, a1', a2' etc. 
  """   
  def get_coefs(self, k, prev_times, h):
    print ("k", k)
    coefs = zeros(6)
    times = (prev_times[0] + h) - array(range(6))*h
    # Fixed step coefficients
    weights = self.get_coefs_dp_dt_at_t0(k, times)
    print ("prev_times", prev_times)
    print ("times", times)
    print weights
    
    coefs[0] = weights[0]
    
    for i in range(1, k+1):
      t = times[i]
      w = weights[i]
      # Coefficients for predictor polynomial 
      print (t,w)
      print shift(self.get_pred_coefs_at_t(k, prev_times, t), 1, cval = 0.0)
      print 
      coefs += w*shift(self.get_pred_coefs_at_t(k, prev_times, t), 1, cval=0.0)

    print coefs
  
  
  # k is degree 
  # times = [t0, t1, ..., t5] is an array of time values. 
  # Suppose that p(t) is the interpolating polynomial passing through
  # (t0, y0) ... (tk, yk). This function evaluates the derivative p'(t) at
  # t = t0
  def get_coefs_dp_dt_at_t0(self, k, times):
    coefs = zeros(6)  
    for i in range(1,k+1):
      coefs += self.get_dd_n(i, times) * self.get_poly_weight_n(i, times)    
    return coefs
    
  # A divided difference [y0, ... , yn] returns a linear comibination
  # of y0, ..., yn. This function returns the weights of that linear combination.
  def get_dd_n(self, n, times):
    n = min(n, 5)
    weights = zeros(6)
    
    for i in range(n+1):
      wi = (times[i] - times[:n+1])
      wi[i] = 1.0
      wi = prod(wi)
      weights[i] = 1.0 / wi
      
    return weights
    
  # Compute the derivative of the polynomial (t - t0)(t - t1) ... (t - t_{n-1}) 
  # and evaluate at t = t0 
  def get_poly_weight_n(self, n, times):
    if n == 1:
      return 1.0
    else :
      return prod(times[0] - times[1:n])
    
      
  
  # Compute the error constant for degree n
  def compute_error_constant(self,n):
    coefs = self.compute_coefficients(n)
    print coefs
    print sum(coefs)
    print coefs[1] + 2.0*coefs[2] + 1.0
    print (0.25)*(coefs[1] + 4.0*coefs[2])
    # Compute
    cn = (1.0 / n) * sum(coefs[1:] * array(range(1,6))**n)
    
    #if n == 1:
    #  cn += 

    
    print array(range(1,7))**n
    print cn
      

bdf_helper = BDF_Helper()
bdf_helper.get_coefs(3, array([6.0, 5.5, 5.0, 4.5, 2.0, 1.0]), 0.5)

"""from pylab import *

ys = array([1.0, 2.0, 0.5, 3.0, 0.0, 0.0])
times = array([5.0, 4.0, 3.0, 2.0, 1.0, 0.0])

ts = linspace(0.0, 6.0, 100)
ps = []
for t in ts:
  ps.append(dot(bdf_helper.get_pred_coefs_at_t(1, times, t), ys))
  
plot(ts, ps)
plot(times, ys, 'ko-')
show()"""
