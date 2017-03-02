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

class BDF_Helper(object):
  
  def __init__(self):
    # Times for last 6 solutions   
    #[t0, t1, t2, t3, t4, t5]
    self.ts = zeros(6)
    

  # Compute time variable coefficients for expression of form 
  # f(t0, y0) = C0*y0 + C1*y1 + ... + CK*yn
  # That is the, coefficients for the time variable nth order BDF method
  def compute_coefficients(self, n):
    weights = zeros(6)
    
    for i in range(1,n+1):
      weights += self.get_dd_n(i) * self.get_poly_weight_n(i)
      
    return weights
    
    
  # The local truncation error LTE is an estimate of how far off the BDF
  # solution is after a time step. For BDF order n it is simply the (n+1)st term
  # of phidot
  def compute_lte_coefficients(self, n):
    return self.get_dd_n(n+1)*self.get_poly_weight_n(n+1)
    
  
  # A divided difference [y0, ... , yn] returns a linear comibination
  # of y0, ..., yn. This function returns the weights of that linear combination.
  def get_dd_n(self, n):
    n = min(n, 5)
    weights = zeros(6)
    
    for i in range(n+1):
      wi = (self.ts[i] - self.ts[:n+1])
      wi[i] = 1.0
      wi = prod(wi)
      weights[i] = 1.0 / wi
      
    return weights
    
  
  # Compute the derivative of the polynomial (t - t0)(t - t1) ... (t - t_{n-1}) 
  # and evaluate at t = t0 (n >= 2)
  def get_poly_weight_n(self, n):
    n = min(n, 6)
    if n == 1:
      return 1.0
    else :
      return prod(self.ts[0] - self.ts[1:n])
      
  
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
bdf_helper.ts = array([6.0, 5.0, 4.0, 3.0, 2.0, 1.0])


bdf_helper.compute_error_constant(2)
#print bdf_helper.compute_coefficients(2)
#print bdf_helper.compute_lte_coefficients(2)
#bdf_helper.ddweights(3)"""
