from dolfin import *
import numpy as np
import matplotlib.pyplot as plt


def Save_XDMF(V, title):
    U = Function(V)
    # input_file = XDMFFile(mesh.mpi_comm(), "solution/RB/%s.xdmf" % title)
    input_file = XDMFFile("solution/%s.xdmf" % title)
    input_file.write(U, "solution")
    input_file.close()


def Load_XDMF(V, title):
    U = Function(V)
    # input_file = XDMFFile(mesh.mpi_comm(), "solution/RB/%s.xdmf" % title)
    input_file = XDMFFile("solution/%s.xdmf" % title)
    input_file.read(U, "solution")
    input_file.close()
    return U


def Save_HDF5(u, mesh, title=None):
    output_file = HDF5File(mesh.mpi_comm(), "solution/%s.h5" % title, "w")
    # input_file = HDF5File("solution/%s.h5" % title, "r")
    output_file.write(u, "solution")
    output_file.close()


def Load_HDF5(V, mesh, title=None):
    U = Function(V)
    input_file = HDF5File(mesh.mpi_comm(), "solution/%s.h5" % title, "r")
    # input_file = HDF5File("solution/%s.h5" % title, "r")
    input_file.read(U, "solution")
    input_file.close()
    return U


""" saving files """


def save_txt(u, title):
    """ saving files """
    coor = mesh.coordinates()
    u_text = []
    u_array = u.vector().get_local()
    print("len %s:%d" % (title, len(u_array)))

    # for i in range(len(u_array)):
    #     u_temp = (coor[i][0], coor[i][1], u_array[i])
    #     print("u(%8g,%8g) = %g" % (coor[i][0], coor[i][1], u_array[i]))
    #     u_text.append(coor[i][0])

    np.savetxt(
        "solution/%s.txt" % title,
        np.array(u11.vector().get_local()),
        #    fmt="%s",
    )
    # np.savetxt("results/periodic_u11_coordinate.txt", np.array(coor))
