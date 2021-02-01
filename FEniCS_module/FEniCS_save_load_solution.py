from dolfin import *
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


'''FEniCS plot '''


def my_FEniCS_plot(number_of_figure, u, title):
    plt.figure(number_of_figure)
    plt.colorbar(plot(u, title='%s' % title))
    # plt.axis('tight')
    # plt.legend()
    # plt.grid(True)
    plt.title(title)
    plt.xlabel('$y_1$')
    plt.ylabel('$y_2$')
    plt.savefig('solution/%s' % title, format='eps')


def my_FEniCS_plot_mode(number_of_figure, u, title, mode):
    plt.figure(number_of_figure)
    plt.colorbar(plot(u, title='%s' % title, mode=mode))
    # plt.axis('tight')
    # plt.legend()
    # plt.grid(True)
    plt.title(title)
    plt.xlabel('$y_1$')
    plt.ylabel('$y_2$')
    plt.savefig('solution/%s' % title, format='eps')