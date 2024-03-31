from Python_module_Quang import *
from ufl import indices
from dolfin import *
from mshr import *
# import math
import meshio
import matplotlib.pyplot as plt
import numpy as np
import time
start = time.time()

filename = "unit_cell.xdmf"
# mesh = load_mesh_RBniCS(filename)


# mesh = meshio.read(
#     filename,  # string, os.PathLike, or a buffer/open file
#     file_format="xdmf",  # optional if filename is a path; inferred from extension
# )

mesh = Mesh()

input_file = HDF5File(mesh.mpi_comm(), "test.h5", "r")
# input_file = HDF5File("solution/%s.h5" % title, "r")
input_file.read(mesh)
input_file.close()


# ------------------------------------------------------------------------------

# mesh_vtk = File("unit_cell.pvd")
# mesh_vtk << mesh
plot(mesh)
plt.show()
