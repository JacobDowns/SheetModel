from dolfin import *
from pylab import *

mesh = Mesh("mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

B = Function(V_cg)
File("B.xml") >> B

coord = mesh.coordinates()
fi = mesh.cells()
v = B.compute_vertex_values(mesh)
vx = coord[:,0]
vy = coord[:,1]

print v.max() - v.min()
print v.min()
print v.max()

#levels = arange(-450,1200,50)
tricontourf(vx, vy, fi, v, 50, cmap=plt.cm.terrain)

colorbar()
show()
 