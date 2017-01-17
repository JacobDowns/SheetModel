import numpy as np

def valley(x,y,bed_para=None): 
	''' bed, ice_thickness=valley(x,y,bed_para)
	A valley shaped topography with fixed surface in which the bed can be changed (with
	bed_para). It is setup such that changing bed_para does not change the outline.
	
	Default bed_para=300/6e3 leads to a glacier which mimics Bench-glacier's bed.
	Decreasing it will lead to thicker ice and eventually an overdeepened bed.
	
	Input:
	x, y coordinates and the bed parameter as defined for the SHMIP exercise, default is 300/6e3

	Output:
	bed and ice thickness in this order

	'''
	# Parameters
	# domain length
	xend = 6.0e3
	# surfaces parameters
	beta = 0.25
	s2 = 100./xend
	min_thick = 1.0
	# bed parameters
	g1 = .5e-6
	alpha = 3.
	bench_para = 300./xend
	if bed_para is None:
		para=bench_para
	else:
		para=bed_para

	# surface
	surf   = 100.*(x+200.)**beta + s2*x - 2.0e10**beta + min_thick
	s_xend = 100.*(xend+200.)**beta + s2*xend - 100.*(200.)**beta + min_thick


	# helper functions
	f_func	= para*x + x**2. * (s_xend-para*6.0e3)/6.0e3**2.	
	f_Bench = bench_para*x + x**2. * (s_xend-bench_para*6.0e3)/6.0e3**2.
	g_func	= 0.5e-6 * abs(y)**3.
	h_func	= (5 - 4.5*x/6.0e3) * (surf-f_func)/(surf-f_Bench)

	#base
	
	bed=f_func+ g_func * h_func
	
	thickness =surf-bed
	thickness[np.where(thickness<0)]=0.0

	return bed, thickness
