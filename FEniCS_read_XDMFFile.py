from dolfin import *
# from Python_module_Quang import FEniCS_plot
# import matplotlib.pyplot as plt
# import numpy

# plt.jet()

# filename_RB = "solution/online_solution_0.xdmf"
# filename_FE = "solution/20210122_2D_elasticity_index_0.h5"

# def test_01():  # working
#     mesh = UnitSquareMesh(10, 10)

#     V = FunctionSpace(mesh, 'P', 1)

#     u = Function(V)
#     u.vector()[:] = 2
#     print(u.vector().get_local())
#     with XDMFFile("test/out.xdmf") as outfile:
#         outfile.write(mesh)
#         outfile.write_checkpoint(u, "u", 0, append=True)

#     mesh2 = Mesh()
#     u = Function(V)
#     with XDMFFile("test/out.xdmf") as infile:
#         infile.read(mesh2)
#         infile.read_checkpoint(u, "u")
#     print(u.vector().get_local())
#     print(mesh.num_cells(), mesh2.num_cells())

# def test_02(filename_RB_checkpoint=None):  # working
#     mesh = Mesh("20210118_data_elastic/elastic_block.xml")
#     V = VectorFunctionSpace(mesh, 'P', 1)
#     # filename_RB = "solution/online_solution_0/solution.xdmf"
#     filename_RB_checkpoint = "solution/online_solution_0/solution_checkpoint.xdmf"

#     # mesh2 = Mesh()
#     u = Function(V)
#     # with XDMFFile(MPI.comm_world, filename_RB_checkpoint) as infile_mesh:
#     #     infile_mesh.read(mesh2)
#     with XDMFFile(MPI.comm_world, filename_RB_checkpoint) as infile_checpoint:
#         infile_checpoint.read_checkpoint(u, "function")
#     # print(u.vector().get_local())
#     plot(u, mode="displacement")
#     # print(mesh.num_cells(), mesh2.num_cells())
#     # f_in = XDMFFile("test/test.xdmf")
#     # f_in.read_checkpoint(f1, "f", 0)
#     return u

# def test_03():  # not working
#     mesh = Mesh("20210118_data_elastic/elastic_block.xml")
#     V = VectorFunctionSpace(mesh, 'P', 1)

#     mesh2 = Mesh()
#     u = Function(V)
#     # with HDF5File(MPI.comm_world, filename_FE) as infile:
#     with HDF5File(mesh.mpi_comm(), filename_FE, "r") as infile:
#         # infile.read(mesh2)
#         infile.read(u, 'solution')
#     LagrangeInterpolator.interpolate(u, mesh)
#     print(u.vector().get_local())
#     print(mesh.num_cells(), mesh2.num_cells())

# def test_04():  # working
#     nx = ny = 8
#     mesh = UnitSquareMesh(nx, ny)
#     V = FunctionSpace(mesh, 'P', 1)

#     alpha = 3
#     beta = 1.2
#     f = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t',
#                    degree=2,
#                    alpha=alpha,
#                    beta=beta,
#                    t=0)
#     f_out = XDMFFile("test/test.xdmf")

#     f_out.write_checkpoint(project(f, V), "f", 0.0, XDMFFile.Encoding.HDF5,
#                            False)

#     for j in range(1, 5):
#         t = float(j / 5)
#         f.t = t

#         f_out.write_checkpoint(project(f, V), "f", t, XDMFFile.Encoding.HDF5,
#                                True)

#     f_out.close()

#     f1 = Function(V)
#     f_in = XDMFFile("test/test.xdmf")

#     f_in.read_checkpoint(f1, "f", 0)
#     f_in.read_checkpoint(f1, "f", 1)
#     f_in.close()


def Read_RB_XDMFFile(V,
                     filename_RB_checkpoint=None,
                     filename_mesh=None,
                     index=None):
    # mesh = Mesh("20210118_data_elastic/elastic_block.xml")
    mesh = Mesh(filename_mesh)
    # V = VectorFunctionSpace(mesh, 'P', 1)
    # filename_RB_checkpoint = "solution/online_solution_0/solution_checkpoint.xdmf"
    # exec(f'u_{index} = k')
    u = Function(V)
    with XDMFFile(MPI.comm_world, filename_RB_checkpoint) as infile_checpoint:
        infile_checpoint.read_checkpoint(u, "function")
        infile_checpoint.close()
    # print(u.vector().get_local())
    # plot(u, mode="displacement")
    # print(mesh.num_cells(), mesh2.num_cells())
    return u


'''
if __name__ == "__main__":
    # test_02()
    # U_0 = Load_HDF5(V=V, mesh=mesh, title='20210122_2D_elasticity_index_0')
    for index in range(5):
        # exec(f'u_{index} = k')
        u = Read_RB_XDMFFile(
            filename_RB_checkpoint=
            "solution/online_solution_%d/solution_checkpoint.xdmf" % index,
            filename_mesh="20210118_data_elastic/elastic_block.xml",
            # index=index
        )
        exec(f'u_{index} = u')
        title = 'u_%s' % index
        plt.figure(index + 1)
        # plot(u, mode="displacement")
        # plt.savefig('solution/%s.png' % title, format='png')
        FEniCS_plot.FEniCS_plot_mode(
            u=u,
            title=title,
            mode="displacement",
            #  number_of_figure=index + 1,
            savefig=True)
    # plt.show(block=False)
'''