"""Unit tests for interpolation using LagrangeInterpolator"""

# Copyright (C) 2014 Mikael Mortensen
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2014-02-18
# Last changed:

import pytest
import numpy
from dolfin import *
import matplotlib.pyplot as plt
plt.jet()
from Python_module_Quang import *


class Quadratic2D(UserExpression):
    def eval(self, values, x):
        values[0] = x[0] * x[0] + x[1] * x[1] + 1.0


class Quadratic3D(UserExpression):
    def eval(self, values, x):
        values[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + 1.0


def test_functional2D():
    """Test integration of function interpolated in non-matching meshes"""
    """LagrangeInterpolator.interpolate(u_project, u_base)"""

    f = Quadratic2D(degree=2)

    # Interpolate quadratic function on course mesh
    mesh0 = UnitSquareMesh(3, 3)
    V0 = FunctionSpace(mesh0, "Lagrange", 2)
    u0 = Function(V0)
    LagrangeInterpolator.interpolate(u0, f)
    plt.figure(1)
    title = "u0_3"
    plot(u0, title=title)

    # Interpolate FE function on finer mesh
    mesh1 = UnitSquareMesh(31, 31)
    V1 = FunctionSpace(mesh1, "Lagrange", 2)
    u1 = Function(V1)
    plt.figure(2)
    title = "u1"
    plot(u1, title=title)
    plt.savefig('%s.png' % title, format='png')

    LagrangeInterpolator.interpolate(u1, u0)
    assert round(assemble(u0 * dx) - assemble(u1 * dx), 10) == 0
    plt.figure(3)
    title = "u0_3"
    plot(u0, wireframe=True, title=title)
    plt.savefig('%s.png' % title, format='png')
    plt.figure(4)
    title = "u1_31"
    plot(u1, title=title)
    plt.savefig('%s.png' % title, format='png')

    # mesh1 = UnitSquareMesh(15, 15)
    # V1 = FunctionSpace(mesh1, "Lagrange", 2)
    # u1 = Function(V1)
    # LagrangeInterpolator.interpolate(u1, u0)
    # assert round(assemble(u0 * dx) - assemble(u1 * dx), 10) == 0
    # plt.figure(4)
    # title = "u0_15"
    # plot(u0, title=title)
    # plt.savefig('%s.png' % title, format='png')


def test_functional3D():
    """Test integration of function interpolated in non-matching meshes"""

    f = Quadratic3D(degree=2)

    # Interpolate quadratic function on course mesh
    mesh0 = UnitCubeMesh(4, 4, 4)
    V0 = FunctionSpace(mesh0, "Lagrange", 2)
    u0 = Function(V0)
    LagrangeInterpolator.interpolate(u0, f)

    # Interpolate FE function on finer mesh
    mesh1 = UnitCubeMesh(11, 11, 11)
    V1 = FunctionSpace(mesh1, "Lagrange", 2)
    u1 = Function(V1)
    LagrangeInterpolator.interpolate(u1, u0)
    assert round(assemble(u0 * dx) - assemble(u1 * dx), 10) == 0

    mesh1 = UnitCubeMesh(10, 11, 10)
    V1 = FunctionSpace(mesh1, "Lagrange", 2)
    u1 = Function(V1)
    LagrangeInterpolator.interpolate(u1, u0)
    assert round(assemble(u0 * dx) - assemble(u1 * dx), 10) == 0


import matplotlib.pyplot as plt

if __name__ == "__main__":
    # import os
    # pytest.main(os.path.abspath(__file__))
    test_functional2D()
    plt.show()