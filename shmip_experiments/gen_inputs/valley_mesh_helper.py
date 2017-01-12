### Valley mesh helper file 

para_bench = 0.05

def surface(x,y):
  return 100.0*(x + 200.0)**(1.0/4.0) + 1.0/60.0*x - 2e10**(1.0/4.0) + 1.0
  
def g(y):
  return 0.5e-6*abs(y)**3
  
def ginv(x):
  return (x/0.5e-6)**(1.0/3.0)
  
def f(x,para):
  return (surface(6e3,0.0) - para*6e3)/6e3**2 * x**2 + para*x
  
def h(x, para):
  return (-4.5*x/6e3 + 5.0) * (surface(x,0.0)-f(x, para)) / (surface(x,0.0)-f(x, para_bench) + 1e-16)

def bed(x,y, para):
  return f(x,para) + g(y)*h(x,para)

def outline(x):
  return ginv( (surface(x,0.0) - f(x,0.05))/(h(x,0.05)+1e-16) )