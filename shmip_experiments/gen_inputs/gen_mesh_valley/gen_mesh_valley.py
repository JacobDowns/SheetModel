from dolfin import *
from mshr import *
import numpy as np
from valley_mesh_helper import *

""" Generate a valley glacier mesh for the SHMIP E experiments. """

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


