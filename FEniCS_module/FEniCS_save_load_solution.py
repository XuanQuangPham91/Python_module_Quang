from dolfin import *
import numpy as np
import meshio

# import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
''' XDMF '''


def save_XDMF(u, mesh, title):
    # u = Function(V)
    # input_file = XDMFFile(mesh.mpi_comm(), "solution/RB/%s.xdmf" % title)
    input_file = XDMFFile(mesh.mpi_comm(), "solution/%s.xdmf" % title)
    # input_file = XDMFFile("solution/%s.xdmf" % title)
    input_file.write_checkpoint(u, "solution")
    input_file.close()


def load_XDMF(V, title):
    u = Function(V)
    # input_file = XDMFFile(mesh.mpi_comm(), "solution/RB/%s.xdmf" % title)
    input_file = XDMFFile("solution/%s.xdmf" % title)
    input_file.read(u, "solution")
    input_file.close()
    return u


def load_mesh_RBniCS(filename):
    """ 
    # Load mesh from xdmf example. 
    * Notice that filename don't have the extension xdmf 
    * Function will create xml file to export mesh as xml
    * Read back and return mesh

    if __name__ == "__main__":
        filename = "solution copy"
        mesh = load_mesh_RBniCS(filename)

        plot(mesh)
        plt.show()
    """
    with meshio.xdmf.TimeSeriesReader(f"{filename}.xdmf") as reader:
        points, cells = reader.read_points_cells()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)

    meshio.Mesh(points, cells).write(f"{filename}.xml")
    mesh = Mesh(f"{filename}.xml")
    return mesh


#------------------------------------------------------------------------------
''' HDF5 '''


def save_HDF5(u, mesh, title=None):
    output_file = HDF5File(mesh.mpi_comm(), "solution/%s.h5" % title, "w")
    # input_file = HDF5File("solution/%s.h5" % title, "r")
    output_file.write(u, "solution")
    output_file.close()


def load_HDF5(V, mesh, title=None):
    u = Function(V)
    input_file = HDF5File(mesh.mpi_comm(), "solution/%s.h5" % title, "r")
    # input_file = HDF5File("solution/%s.h5" % title, "r")
    input_file.read(u, "solution")
    input_file.close()
    return u


#------------------------------------------------------------------------------
''' vtk '''


def save_vtk(u, title):
    # Save solution to file in VTK format
    File('solution/%s.pvd' % title) << u


#------------------------------------------------------------------------------
''' txt '''


def save_u_txt(u, title):
    u_array = u.vector().get_local()
    print("len of %s: %d" % (title, len(u_array)))
    np.savetxt(
        "solution/%s.txt" % title,
        np.array(u_array),
        #    fmt="%s",
    )

    # u_text = []
    # coor = mesh.coordinates()
    # for i in range(len(u_array)):
    #     u_temp = (coor[i][0], coor[i][1], u_array[i])
    #     print("u(%8g,%8g) = %g" % (coor[i][0], coor[i][1], u_array[i]))
    #     u_text.append(coor[i][0])
