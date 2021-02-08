""" FEniCS auxiliary codes """

from dolfin import *
# from rbnics import *
# from Python_module_Quang import *
# import clear_offline_data
# import matplotlib.pyplot as plt
import numpy as np
# import math


def normalize_solution(u):
    "Normalize u: return u divided by max(u)"
    u_array = u.vector().get_local()
    u_max = np.max(np.abs(u_array))
    u_array /= u_max
    u.vector()[:] = u_array
    # u.vector().set_local(u_array)  # alternative
    return u


def absolute_error(u_e, u, space="VectorFunctionSpace"):
    V = u.function_space()
    mesh = V.mesh()
    deggre = V.ufl_element().degree()
    if space == "VectorFunctionSpace":
        W = VectorFunctionSpace(mesh, 'P', degree=3)
    elif space == "FunctionSpace":
        W = FunctionSpace(mesh, 'P', 1)

    # u_e_W = LagrangeInterpolator.interpolate(u_e, W)
    # u_W = LagrangeInterpolator.interpolate(u, W)
    u_e_W = interpolate(u_e, W)
    u_W = interpolate(u, W)
    # u_e_W = u_e
    # u_W = u
    e_W = Function(W)
    e_W.vector()[:] = np.absolute(u_e_W.vector().get_local() -
                                  u_W.vector().get_local())
    # e_W.vector()[:] = np.absolute(u_e_W.vector().get_local() - u_W.vector().get_local())/u_e_W.vector().get_local()
    return e_W


def relative_error(u_e, u, space="VectorFunctionSpace"):
    V = u.function_space()
    mesh = V.mesh()
    deggre = V.ufl_element().degree()
    if space == "VectorFunctionSpace":
        W = VectorFunctionSpace(mesh, 'P', degree=3)
    elif space == "FunctionSpace":
        W = FunctionSpace(mesh, 'P', 3)
    # u_e_W = LagrangeInterpolator.interpolate(u_e, W)
    # u_W = LagrangeInterpolator.interpolate(u, W)
    u_e_W = interpolate(u_e, W)
    u_W = interpolate(u, W)
    e_W = Function(W)
    # e_W.vector()[:] = np.absolute(u_e_W.vector().get_local() - u_W.vector().get_local())
    e_W.vector()[:] = np.absolute(1 - (u_W.vector().get_local() /
                                       u_e_W.vector().get_local()))
    return e_W