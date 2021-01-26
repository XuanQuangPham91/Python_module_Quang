from dolfin import *


def Save_XDMF(V, title):
    U = Function(V)
    # input_file = XDMFFile(mesh.mpi_comm(), "solution/RB/%s.xdmf" % title)
    input_file = XDMFFile("solution/RB/%s.xdmf" % title)
    input_file.write(U, "solution")
    input_file.close()


def Load_XDMF(V, title):
    U = Function(V)
    # input_file = XDMFFile(mesh.mpi_comm(), "solution/RB/%s.xdmf" % title)
    input_file = XDMFFile("solution/RB/%s.xdmf" % title)
    input_file.read(U, "solution")
    input_file.close()
    return U


def Save_HDF5(V, mesh, title=None):
    # Load solution
    U = Function(V)
    input_file = HDF5File(mesh.mpi_comm(), "solution/RB/%s.h5" % title, "r")
    # input_file = HDF5File("solution/%s.h5" % title, "r")
    input_file.read(U, "solution")
    input_file.close()


def Load_HDF5(V, mesh, title=None):
    # Load solution
    U = Function(V)
    input_file = HDF5File(mesh.mpi_comm(), "solution/RB/%s.h5" % title, "r")
    # input_file = HDF5File("solution/%s.h5" % title, "r")
    input_file.read(U, "solution")
    input_file.close()
    return U