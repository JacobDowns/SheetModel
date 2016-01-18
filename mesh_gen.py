# -*- coding: utf-8 -*-
"""
Create a rectangular test mesh.
"""

from varglas.meshing           import MeshGenerator
from varglas.io                import DataInput
from varglas.data.data_factory import DataFactory
import numpy as np


b = DataFactory.get_bamber()
db = DataInput(b, gen_space = False) 
m = MeshGenerator(db, 'mesh', 'meshes/')

xs = [0.0, 60e3, 60e3, 0.0]
ys = [0.0, 0.0, 20e3, 20e3] 

domain_coordinates = np.transpose([xs, ys])

m.set_contour(domain_coordinates)
m.write_gmsh_contour(500)
m.close_file()
