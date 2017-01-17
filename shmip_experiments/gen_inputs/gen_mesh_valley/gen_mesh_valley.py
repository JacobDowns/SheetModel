#from dolfin import *
#from mshr import *
#import numpy as np
from valley_outline import *
from pylab import *
from mshr import *

""" Generate a valley glacier mesh for the SHMIP E experiments. """


xs, ys = valley_outline()

print xs 
print ys

print len(xs)
plot(xs[:1000], ys[:1000], 'ko-', ms = 2)
#show()

savefig('out.png')



quit()

# Domains length
L = 6000.0
# Halfwidth of the domain
xs = np.linspace(0.0, L, 300)
ys = -outline(xs)


# Build the polygonal domain
points = []

for i in range(len(xs)):
  points.append(Point(xs[i], ys[i]))

xs = xs[::-1][1:]
ys = ys[::-1][1:]

for i in range(len(xs)):
  points.append(Point(xs[i], -ys[i]))
  
# Create the mesh
domain = Polygon(points)
mesh = generate_mesh(domain, 200) 
plot(mesh, interactive = True)
File('mesh_valley.xml') << mesh


