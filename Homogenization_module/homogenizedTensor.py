# %%
""" 2D Linear elastic with isotropic material """
""" Note: 2020.11.23
* Test and reference to
  * Macedo2018. Configuration: 1 fiber in the middle and 4 corners
  * Oliveira2009. Configuration: 1 fiber in 4 corners
* Each cofiguration was prepared under comment form. Can be used by decomment.
* Pictures will be saved under eps type.
* Saving Paraview file (.pvd) was prepared under comments and be considered as option.
"""

# import libraries
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
import time

# import self-packages
from Python_module_Quang import *

# %%
# Define initial paramemters --------------------------------------------------
# V_f = 0.47  # volume of fiber - J.A Oliveira

# %%
""" Numerical integration of Homogenized elastic tensor (efficient tensor) """
""" Compute gradient of u_kl """

class homogenizedTensor():

    # Default initialization of members
    def __init__(self, material_properties, u11, u22, u12, **kwargs):
        # Call the standard initialization
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        assert "mesh" in kwargs
        self.subdomains, self.boundaries, self.mesh = (
            kwargs["subdomains"],
            kwargs["boundaries"],
            kwargs["mesh"],
        )
        # self.u = TrialFunction(V)
        # self.v = TestFunction(V)
        self.dx = Measure("dx")(subdomain_data=self.subdomains)
        self.ds = Measure("ds")(subdomain_data=self.boundaries)

        # initialize the cell solutions
        self.u11 = u11
        self.u22 = u22
        self.u12 = u12

        # set E & nu
        self.material_properties = material_properties
        self.E_f = material_properties[0]  # MPa
        self.nu_f = material_properties[1]
        self.E_m = material_properties[2]  # MPa
        self.nu_m = material_properties[3]
        # self.E = [self.E_m, self.E_f]
        # self.nu = [self.nu_m, self.nu_f]

        # Lame's constants
        # self.lambda_1 = [
        #     self.E_f * self.nu_f / ((1.0 + self.nu_f) *
        #                             (1.0 - 2.0 * self.nu_f)),
        #     self.E_m * self.nu_m / ((1.0 + self.nu_m) *
        #                             (1.0 - 2.0 * self.nu_m)),
        # ]
        # self.lambda_2 = [
        #     self.E_f / (2.0 * (1.0 + self.nu_f)), \
        #     self.E_m / (2.0 * (1.0 + self.nu_m))
        #     ]

        self.model = "plane_strain"
        if self.model == "plane_strain":
            # define C_ij tensor
            self.C1111 = self.C2222 = [
                self.E_f * (1 - self.nu_f) / ((1.0 + self.nu_f) *(1.0 - 2.0 * self.nu_f)), \
                self.E_m * (1 - self.nu_m) / ((1.0 + self.nu_m) * (1.0 - 2.0 * self.nu_m))
            ]
            self.C1122 = self.C2211 = [
                self.E_f * self.nu_f / ((1.0 + self.nu_f) * (1.0 - 2.0 * self.nu_f)),
                self.E_m * self.nu_m / ((1.0 + self.nu_m) * (1.0 - 2.0 * self.nu_m)),
            ]
            self.C1212 = [
                self.E_f / (2 * (1.0 + self.nu_f)),
                self.E_m / (2 * (1.0 + self.nu_m))
            ]
        elif self.model == "plane_stress":
            self.C1111 = self.C2222 = [
                self.E_f / (1.0 - np.power(self.nu_f, 2)), \
                self.E_m / (1.0 - np.power(self.nu_m, 2))
            ]
            self.C1122 = self.C2211 = [
                self.E_f * self.nu_f / (1.0 - np.power(self.nu_f, 2)),
                self.E_m * self.nu_m / (1.0 - np.power(self.nu_m, 2)),
            ]
            self.C1212 = [
                self.E_f / (2 * (1.0 + self.nu_f)),
                self.E_m / (2 * (1.0 + self.nu_m))
            ]

        # C_0 = as_matrix([
        #     [self.C1111[0], self.C1122[0], 0.],  #
        #     [self.C2211[0], self.C2222[0], 0.],  #
        #     [0., 0., self.C1212[0]]
        # ])
        # C_1 = as_matrix([
        #     [self.C1111[1], self.C1122[1], 0.],  #
        #     [self.C2211[1], self.C2222[1], 0.],  #
        #     [0., 0., self.C1212[1]]
        # ])

        # call homogenizedTensor problem
        self.A_ij, self.E_hom, self.nu_hom = self._homogenizedTensor()

    def _homogenizedTensor(self):
        # Default initialization of members
        dx = self.dx
        ds = self.ds
        u11, u22, u12 = self.u11, self.u22, self.u12
        C1111, C2222 = self.C1111, self.C2222
        C1122, C2211 = self.C1122, self.C2211
        C1212 = self.C1212
        msh = self.mesh
        # Create TensorFunctionSpace for grad of u11, u22, u12
        V_g = TensorFunctionSpace(msh,
                                  "Lagrange",
                                  1,
                                  constrained_domain=PeriodicBoundary())
        V = FunctionSpace(msh,
                          "Lagrange",
                          1,
                          constrained_domain=PeriodicBoundary())

        # Calculate gradient of u11, u22, u12 ----------------------------------
        # for w11 (or u11)
        time_grad_u11_s = time.time()
        grad_u11 = project(grad(u11), V_g)
        grad_u11_1_y1, grad_u11_1_y2, grad_u11_2_y1, grad_u11_2_y2 = \
            grad_u11.split(deepcopy=True)  # extract
        time_grad_u11_e = time.time()

        # for w22 (or u22)
        time_grad_u22_s = time.time()
        grad_u22 = project(grad(u22), V_g)
        grad_u22_1_y1, grad_u22_1_y2, grad_u22_2_y1, grad_u22_2_y2 = \
            grad_u22.split(deepcopy=True)
        time_grad_u22_e = time.time()

        # for w12 (or u12)
        time_grad_u12_s = time.time()
        grad_u12 = project(grad(u12), V_g)
        grad_u12_1_y1, grad_u12_1_y2, grad_u12_2_y1, grad_u12_2_y2 = \
            grad_u12.split(deepcopy=True)
        time_grad_u12_e = time.time()

        # ----------------------------------------------------------------------
        """ Subtitute into equations of homogenized tensor components """

        ## k = l = 1 -----------------------------------------------------------
        time_A1111_s = time.time()
        A1111_projected_1 = project(
            C1111[0] - C1111[0] * grad_u11_1_y1 - C1122[0] * grad_u11_2_y2, V=V)
        A1111_projected_2 = project(
            C1111[1] - C1111[1] * grad_u11_1_y1 - C1122[1] * grad_u11_2_y2, V=V)
        A1111_assembled = assemble(A1111_projected_1 * dx(1) +
                                   A1111_projected_2 * dx(2) +
                                   A1111_projected_2 * dx(3) +
                                   A1111_projected_2 * dx(4) +
                                   A1111_projected_2 * dx(5) +
                                   A1111_projected_2 * dx(6) +
                                   A1111_projected_2 * dx(7) +
                                   A1111_projected_2 * dx(8) +
                                   A1111_projected_2 * dx(9) +
                                   A1111_projected_2 * dx(10))
        time_A1111_e = time.time()

        time_A2211_s = time.time()
        A2211_projected_1 = project(
            C2211[0] - C2211[0] * grad_u11_1_y1 - C2222[0] * grad_u11_2_y2, V=V)
        A2211_projected_2 = project(
            C2211[1] - C2211[1] * grad_u11_1_y1 - C2222[1] * grad_u11_2_y2, V=V)
        A2211_assembled = assemble(A2211_projected_1 * dx(1) +
                                   A2211_projected_2 * dx(2) +
                                   A2211_projected_2 * dx(3) +
                                   A2211_projected_2 * dx(4) +
                                   A2211_projected_2 * dx(5) +
                                   A2211_projected_2 * dx(6) +
                                   A2211_projected_2 * dx(7) +
                                   A2211_projected_2 * dx(8) +
                                   A2211_projected_2 * dx(9) +
                                   A2211_projected_2 * dx(10))
        time_A2211_e = time.time()

        ## k = l = 2 -----------------------------------------------------------
        time_A1122_s = time.time()
        A1122_projected_1 = project(
            C1122[0] - C1111[0] * grad_u22_1_y1 - C1122[0] * grad_u22_2_y2, V=V)
        A1122_projected_2 = project(
            C1122[1] - C1111[1] * grad_u22_1_y1 - C1122[1] * grad_u22_2_y2, V=V)
        A1122_assembled = assemble(A1122_projected_1 * dx(1) +
                                   A1122_projected_2 * dx(2) +
                                   A1122_projected_2 * dx(3) +
                                   A1122_projected_2 * dx(4) +
                                   A1122_projected_2 * dx(5) +
                                   A1122_projected_2 * dx(6) +
                                   A1122_projected_2 * dx(7) +
                                   A1122_projected_2 * dx(8) +
                                   A1122_projected_2 * dx(9) +
                                   A1122_projected_2 * dx(10))
        time_A1122_e = time.time()

        time_A2222_s = time.time()
        A2222_projected_1 = project(
            C2222[0] - C2211[0] * grad_u22_1_y1 - C2222[0] * grad_u22_2_y2, V=V)
        A2222_projected_2 = project(
            C2222[1] - C2211[1] * grad_u22_1_y1 - C2222[1] * grad_u22_2_y2, V=V)
        A2222_assembled = assemble(A2222_projected_1 * dx(1) +
                                   A2222_projected_2 * dx(2) +
                                   A2222_projected_2 * dx(3) +
                                   A2222_projected_2 * dx(4) +
                                   A2222_projected_2 * dx(5) +
                                   A2222_projected_2 * dx(6) +
                                   A2222_projected_2 * dx(7) +
                                   A2222_projected_2 * dx(8) +
                                   A2222_projected_2 * dx(9) +
                                   A2222_projected_2 * dx(10))
        time_A2222_e = time.time()

        # k =1, l = 2 ----------------------------------------------------------
        time_A1212_s = time.time()
        A1212_projected_1 = project(
            C1212[0] * (1 - grad_u12_1_y2 - grad_u12_2_y1), V=V)
        A1212_projected_2 = project(
            C1212[1] * (1 - grad_u12_1_y2 - grad_u12_2_y1), V=V)
        A1212_assembled = assemble(A1212_projected_1 * dx(1) +
                                   A1212_projected_2 * dx(2) +
                                   A1212_projected_2 * dx(3) +
                                   A1212_projected_2 * dx(4) +
                                   A1212_projected_2 * dx(5) +
                                   A1212_projected_2 * dx(6) +
                                   A1212_projected_2 * dx(7) +
                                   A1212_projected_2 * dx(8) +
                                   A1212_projected_2 * dx(9) +
                                   A1212_projected_2 * dx(10))
        # A1212_assembled = assemble(A1212_projected_1 * dx(1))\
        #                 + assemble(A1212_projected_2 * dx(2))\
        #                 + assemble(A1212_projected_2 * dx(3))\
        #                 + assemble(A1212_projected_2 * dx(4))\
        #                 + assemble(A1212_projected_2 * dx(5))\
        #                 + assemble(A1212_projected_2 * dx(6))\
        #                 + assemble(A1212_projected_2 * dx(7))\
        #                 + assemble(A1212_projected_2 * dx(8))\
        #                 + assemble(A1212_projected_2 * dx(9))\
        #                 + assemble(A1212_projected_2 * dx(10))
        # a1 = [assemble(A1212_projected_1 * dx(i+1)) for i in range(1)]
        # a2 = [assemble(A1212_projected_2 * dx(i+1)) for i in range(1,10)]
        # print(f"assemble1: {a1}")
        # print(f"assemble2: {a2}")
        # print(f"assemble2 total: {np.sum(a2)}")
        # print(f"assemble total: {a1 + np.sum(a2)}")
        time_A1212_e = time.time()
        # print('The homogenized coefficient A1212: ', A1212_assembled / 1e9)

        #=======================================================================
        #=======================================================================
        print('Computational time - calculate gradient of u11, u22, u12: ')
        self.time_total_grad = (time_grad_u12_e - time_grad_u11_s)
        print("""
            time_grad_u11: {} \n
            time_grad_u22: {} \n
            time_grad_u12: {} \n
            time_grad_total: {} \n
        """.format(time_grad_u11_e - time_grad_u11_s,
                   time_grad_u22_e - time_grad_u22_s,
                   time_grad_u12_e - time_grad_u12_s,
                   self.time_total_grad))
        #=======================================================================
        print('Computational time - homogenized coefficient: ')
        self.time_total_assemble = (time_A1212_e - time_A1111_s)
        print("""
            time_A1111: {} \n
            time_A2211: {} \n
            time_A1122: {} \n
            time_A2222: {} \n
            time_A1212: {} \n
            time_total_assemble: {} \n
        """.format(time_A1111_e - time_A1111_s, time_A2211_e - time_A2211_s,
                   time_A1122_e - time_A1122_s, time_A2222_e - time_A2222_s,
                   time_A1212_e - time_A1212_s, self.time_total_assemble))

        #=======================================================================
        # The homogenized tensor
        A_ij = np.array([
            [A1111_assembled, A1122_assembled, 0],
            [A2211_assembled, A2222_assembled, 0],
            [0, 0, A1212_assembled],
        ])
        print("The homogenized effective coefficient tensor (GPa): \n",
              A_ij / 1e3)

        lambda_hom = A_ij[0, 1]
        mu_hom = A_ij[2, 2]

        # ----------------------------------------------------------------------
        """ The homogenized material properties"""
        if self.model == "plane_strain":
            # ...
            # ==================================================================
            # ==================================================================
            # E_hom = mu_hom * (3 * lambda_hom + 2 * mu_hom) / (lambda_hom + mu_hom)
            # nu_hom = lambda_hom / (2 * (lambda_hom + 2*mu_hom))
            E_hom = A_ij[0, 0] - 2*(A_ij[0, 1])/(A_ij[0, 0]+A_ij[0, 1])
            nu_hom = (A_ij[0, 1])/(A_ij[0, 0]+A_ij[0, 1])
            # References
            # [1] Penta, R. and Gerisch, A., 2017. The asymptotic homogenization 
            # elasticity tensor properties for composites with material discontinuities.
            # Continuum Mechanics and Thermodynamics, 29, pp.187-206.
            # ==================================================================
            # Using this one, do not delete ====================================
            # nu_hom = A1122_assembled / (A1122_assembled + A1111_assembled)
            # E_hom = 2 * A1212_assembled * (1 + nu_hom)
            # ==================================================================
            # ==================================================================
        elif self.model == "plane_stress":
            #     # plain stress: eq reference to JinheeLee1996
            nu_hom = A1122_assembled / A1111_assembled
            E_hom = A1111_assembled * (1 - np.power(nu_hom, 2))
            # E_hom = A1212_assembled * 2 * (1 + nu_hom)

        # self.inv_A_ij = np.linalg.inv(A_ij)
        # E_hom = 1 / self.inv_A_ij[0, 0]
        # nu_hom = -self.inv_A_ij[1, 0] * E_hom

        # E_hom = np.round(E_hom / 1e3, 2)
        # nu_hom = np.round(nu_hom, 2)
        # print(f"E_hom: {E_hom} (GPa)")
        # print(f"nu_hom: {nu_hom}")

        print(f"E_hom: {E_hom} (MPa)")
        print(f"nu_hom: {nu_hom}")
        return A_ij, E_hom, nu_hom


class homogenizedTensor3D():
    """
    Homogenized tensor for 3D without pre-computed gradients of correctors
        - Linear elasticity or vector-valued problem
    """
    # Default initialization of members
    def __init__(self, material_properties, u11, u22, u33, u12, u13, u23, n_subs=10, **kwargs):
        """
        n_subs: number of subdomains
        """
        # Call the standard initialization
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        assert "mesh" in kwargs
        self.subdomains, self.boundaries, self.mesh = (
            kwargs["subdomains"],
            kwargs["boundaries"],
            kwargs["mesh"],
        )
        # self.u = TrialFunction(V)
        # self.v = TestFunction(V)
        self.dx = Measure("dx")(domain=self.mesh, subdomain_data=self.subdomains)
        self.ds = Measure("ds")(domain=self.mesh, subdomain_data=self.boundaries)
        self.n_subs = n_subs

        # initialize the cell solutions
        self.u11, self.u22, self.u33 = u11, u22, u33
        self.u12, self.u13, self.u23 = u12, u13, u23

        # set E & nu
        self.material_properties = material_properties
        self.E_f = material_properties[0]  # MPa
        self.nu_f = material_properties[1]
        self.E_m = material_properties[2]  # MPa
        self.nu_m = material_properties[3]
        # self.E = [self.E_m, self.E_f]
        # self.nu = [self.nu_m, self.nu_f]

        # Lame's constants
        self.lambda_1 = [
            self.E_f * self.nu_f / ((1.0 + self.nu_f) * (1.0 - 2.0 * self.nu_f)),
            self.E_m * self.nu_m / ((1.0 + self.nu_m) * (1.0 - 2.0 * self.nu_m)),
        ]
        self.lambda_2 = [self.E_f / (2.0 * (1.0 + self.nu_f)),
                         self.E_m / (2.0 * (1.0 + self.nu_m)),]
        self.C1111 = [self.lambda_1[0] + 2 * self.lambda_2[0],
                      self.lambda_1[1] + 2 * self.lambda_2[1],]
        self.C2222 = self.C3333 = self.C1111
        self.C1122 = [self.lambda_1[0],
                      self.lambda_1[1],]
        self.C2211 = self.C2233 = self.C1133 = self.C3311 = self.C3322 = self.C1122
        self.C2323 = [self.lambda_2[0],
                      self.lambda_2[1],]
        self.C1212 = self.C2323
        self.C1313 = self.C2323
        # call homogenizedTensor problem
        self.A_ij, self.material_properties_homogenized = self._homogenizedTensor()

    def _homogenizedTensor(self):
        # Default initialization of members
        dx = self.dx
        ds = self.ds
        u11, u22, u33 = self.u11, self.u22, self.u33
        u12, u13, u23 = self.u12, self.u13, self.u23
        C1111, C2222, C3333 = self.C1111, self.C2222, self.C3333
        C1122, C2211 = self.C1122, self.C2211
        C1133, C3311 = self.C1133, self.C3311
        C1212, C1313, C2323 = self.C1212, self.C1313, self.C2323
        C2233, C3322 = self.C2233, self.C3322
        msh = self.mesh
        # Create TensorFunctionSpace for grad of u11, u22, u12
        self.V_g = TensorFunctionSpace(msh,
                                  "Lagrange",
                                  1,
                                  constrained_domain=PeriodicBoundary3D())
        self.V = FunctionSpace(msh,
                          "Lagrange",
                          1,
                          constrained_domain=PeriodicBoundary3D())

        # Calculate gradient of u11, u22, u12 ----------------------------------
        # for w11 (or u11)
        time_grad_u11_s = time.time()
        grad_u11 = project(grad(u11), self.V_g)
        grad_u11_1_y1, grad_u11_1_y2, grad_u11_1_y3, \
            grad_u11_2_y1, grad_u11_2_y2, grad_u11_2_y3, \
            grad_u11_3_y1, grad_u11_3_y2, grad_u11_3_y3 = grad_u11.split(deepcopy=True) # extract
        time_grad_u11_e = time.time()

        # for w22 (or u22)
        time_grad_u22_s = time.time()
        grad_u22 = project(grad(u22), self.V_g)
        grad_u22_1_y1, grad_u22_1_y2, grad_u22_1_y3,\
            grad_u22_2_y1, grad_u22_2_y2, grad_u22_2_y3,\
            grad_u22_3_y1, grad_u22_3_y2, grad_u22_3_y3 = grad_u22.split(deepcopy=True)
        time_grad_u22_e = time.time()

        # for w33 (or u12)
        time_grad_u33_s = time.time()
        grad_u33 = project(grad(u33), self.V_g)
        grad_u33_1_y1, grad_u33_1_y2, grad_u33_1_y3,\
            grad_u33_2_y1, grad_u33_2_y2, grad_u33_2_y3,\
            grad_u33_3_y1, grad_u33_3_y2, grad_u33_3_y3 = grad_u33.split(deepcopy=True)
        time_grad_u33_e = time.time()

        # for w23 (or u23)
        time_grad_u23_s = time.time()
        grad_u23 = project(grad(u23), self.V_g)
        grad_u23_1_y1, grad_u23_1_y2, grad_u23_1_y3,\
            grad_u23_2_y1, grad_u23_2_y2, grad_u23_2_y3,\
            grad_u23_3_y1, grad_u23_3_y2, grad_u23_3_y3 = grad_u23.split(deepcopy=True)
        time_grad_u23_e = time.time()

        # for w13 (or u13)
        time_grad_u13_s = time.time()
        grad_u13 = project(grad(u13), self.V_g)
        grad_u13_1_y1, grad_u13_1_y2, grad_u13_1_y3,\
            grad_u13_2_y1, grad_u13_2_y2, grad_u13_2_y3,\
            grad_u13_3_y1, grad_u13_3_y2, grad_u13_3_y3 = grad_u13.split(deepcopy=True)
        time_grad_u13_e = time.time()

        # for w12 (or u12)
        time_grad_u12_s = time.time()
        grad_u12 = project(grad(u12), self.V_g)
        grad_u12_1_y1, grad_u12_1_y2, grad_u12_1_y3,\
            grad_u12_2_y1, grad_u12_2_y2, grad_u12_2_y3,\
            grad_u12_3_y1, grad_u12_3_y2, grad_u12_3_y3 = grad_u12.split(deepcopy=True)
        time_grad_u12_e = time.time()

        #=======================================================================
        """ Subtitute into equations of homogenized tensor components """
        #=======================================================================
        """group 1: k = l = 1"""
        time_A1111_s = time.time()
        A1111_projected_1 = project(C1111[0]
                                    - C1111[0] * grad_u11_1_y1
                                    - C1122[0] * grad_u11_2_y2
                                    - C1133[0] * grad_u11_3_y3, V=self.V)
        A1111_projected_2 = project(C1111[1] 
                                    - C1111[1] * grad_u11_1_y1 
                                    - C1122[1] * grad_u11_2_y2
                                    - C1133[1] * grad_u11_3_y3, V=self.V)
        # A1111_assembled = assemble(
        #     A1111_projected_1 * dx(1) +
        #     A1111_projected_2 * dx(2) +
        #     A1111_projected_2 * dx(3) +
        #     A1111_projected_2 * dx(4) +
        #     A1111_projected_2 * dx(5) +
        #     A1111_projected_2 * dx(6) +
        #     A1111_projected_2 * dx(7) +
        #     A1111_projected_2 * dx(8) +
        #     A1111_projected_2 * dx(9) +
        #     A1111_projected_2 * dx(10))
        A1111_assembled = assemble(A1111_projected_1 * dx(1) +
                                   sum([A1111_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        # A1111_projected = [project(
        #     C1111[i] 
        #     - C1111[i] * grad_u11_1_y1 
        #     - C1122[i] * grad_u11_2_y2
        #     - C1133[i] * grad_u11_3_y3, V=self.V) for i in range(2)]
        # A1111_assembled = assemble(A1111_projected[0] * dx(1) +
        #                            sum([A1111_projected[1] * dx(i+1) for i in range(1, self.n_subs)]))
        time_A1111_e = time.time()

        time_A2211_s = time.time()
        A2211_projected_1 = project(C2211[0] 
                                    - C2211[0] * grad_u11_1_y1
                                    - C2222[0] * grad_u11_2_y2
                                    - C2233[0] * grad_u11_3_y3, V=self.V)
        A2211_projected_2 = project(C2211[1]
                                    - C2211[1] * grad_u11_1_y1
                                    - C2222[1] * grad_u11_2_y2
                                    - C2233[1] * grad_u11_3_y3, V=self.V)
        A2211_assembled = assemble(A2211_projected_1 * dx(1) +
                                   sum([A2211_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A2211_e = time.time()

        time_A3311_s = time.time()
        A3311_projected_1 = project(C3311[0]
                                    - C3311[0] * grad_u11_1_y1
                                    - C3322[0] * grad_u11_2_y2
                                    - C3333[0] * grad_u11_3_y3, V=self.V)
        A3311_projected_2 = project(C3311[1]
                                    - C3311[1] * grad_u11_1_y1
                                    - C3322[1] * grad_u11_2_y2
                                    - C3333[1] * grad_u11_3_y3, V=self.V)
        A3311_assembled = assemble(A3311_projected_1 * dx(1) +
                                   sum([A3311_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A3311_e = time.time()

        #=======================================================================
        """ group 2: k = l = 2 """
        time_A1122_s = time.time()
        A1122_projected_1 = project(
            C1122[0] - C1111[0] * grad_u22_1_y1 - C1122[0] * grad_u22_2_y2, V=self.V)
        A1122_projected_2 = project(
            C1122[1] - C1111[1] * grad_u22_1_y1 - C1122[1] * grad_u22_2_y2, V=self.V)
        A1122_assembled = assemble(A1122_projected_1 * dx(1) +
                                   sum([A1122_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A1122_e = time.time()

        time_A2222_s = time.time()
        A2222_projected_1 = project(
            C2222[0] - C2211[0] * grad_u22_1_y1 - C2222[0] * grad_u22_2_y2, V=self.V)
        A2222_projected_2 = project(
            C2222[1] - C2211[1] * grad_u22_1_y1 - C2222[1] * grad_u22_2_y2, V=self.V)
        A2222_assembled = assemble(A2222_projected_1 * dx(1) +
                                   sum([A2222_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A2222_e = time.time()

        time_A3322_s = time.time()
        A3322_projected_1 = project(C3322[0]
                                - C3311[0] * grad_u22_1_y1
                                - C3322[0] * grad_u22_2_y2
                                - C3333[0] * grad_u22_3_y3, V=self.V)
        A3322_projected_2 = project(C3322[1]
                                - C3311[1] * grad_u22_1_y1
                                - C3322[1] * grad_u22_2_y2
                                - C3333[1] * grad_u22_3_y3, V=self.V)
        A3322_assembled = assemble(A3322_projected_1 * dx(1) +
                                   sum([A3322_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A3322_e = time.time()

        #=======================================================================
        """ Group 3: k = l = 3 """
        time_A1133_s = time.time()
        A1133_projected_1 = project(
            C1122[0] - C1111[0] * grad_u22_1_y1 - C1122[0] * grad_u22_2_y2, V=self.V)
        A1133_projected_2 = project(
            C1122[1] - C1111[1] * grad_u22_1_y1 - C1122[1] * grad_u22_2_y2, V=self.V)
        A1133_assembled = assemble(A1133_projected_1 * dx(1) +
                                   sum([A1133_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A1133_e = time.time()

        time_A2233_s = time.time()
        A2233_projected_1 = project(C2233[0]
                                - C2211[0] * grad_u33_1_y1
                                - C2222[0] * grad_u33_2_y2
                                - C2233[0] * grad_u33_3_y3, V=self.V)
        A2233_projected_2 = project(C2233[1]
                                - C2211[1] * grad_u33_1_y1
                                - C2222[1] * grad_u33_2_y2
                                - C2233[1] * grad_u33_3_y3, V=self.V)
        A2233_assembled = assemble(A2233_projected_1 * dx(1) +
                                   sum([A2233_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A2233_e = time.time()

        time_A3333_s = time.time()
        A3333_projected_1 = project(C3333[0]
                                - C3311[0] * grad_u33_1_y1
                                - C3322[0] * grad_u33_2_y2
                                - C3333[0] * grad_u33_3_y3, V=self.V)
        A3333_projected_2 = project(C3333[1]
                                - C3311[1] * grad_u33_1_y1
                                - C3322[1] * grad_u33_2_y2
                                - C3333[1] * grad_u33_3_y3, V=self.V)
        A3333_assembled = assemble(A3333_projected_1 * dx(1) +
                                   sum([A3333_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A3333_e = time.time()

        #=======================================================================
        """ Group 4:
        # k =2, l = 3 
        # k =1, l = 3
        # k =1, l = 2
        """ 
        time_A2323_s = time.time()
        A2323_projected_1 = project(C2323[0] * (1 - grad_u23_2_y3 - grad_u23_3_y2), V=self.V)
        A2323_projected_2 = project(C2323[1] * (1 - grad_u23_2_y3 - grad_u23_3_y2), V=self.V)
        A2323_assembled = assemble(A2323_projected_1 * dx(1) +
                                   sum([A2323_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A2323_e = time.time()

        time_A1313_s = time.time()
        A1313_projected_1 = project(C1313[0] * (1 - grad_u13_1_y3 - grad_u13_3_y1),V=self.V)
        A1313_projected_2 = project(C1313[1] * (1 - grad_u13_1_y3 - grad_u13_3_y1),V=self.V)
        A1313_assembled = assemble(A1313_projected_1 * dx(1) +
                                   sum([A1313_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A1313_e = time.time()

        time_A1212_s = time.time()
        A1212_projected_1 = project(C1212[0] * (1 - grad_u12_1_y2 - grad_u12_2_y1), V=self.V)
        A1212_projected_2 = project(C1212[1] * (1 - grad_u12_1_y2 - grad_u12_2_y1), V=self.V)
        A1212_assembled = assemble(A1212_projected_1 * dx(1) +
                                   sum([A1212_projected_2 * dx(i+1) for i in range(1, self.n_subs)]))
        time_A1212_e = time.time()
        # print('The homogenized coefficient A1212: ', A1212_assembled / 1e9)

        #=======================================================================
        #=======================================================================
        print('Computational time - calculate gradient of u11, u22, u12: ')
        self.time_total_grad = (time_grad_u12_e - time_grad_u11_s)
        print("""
            time_grad_u11: {} \n
            time_grad_u22: {} \n
            time_grad_u33: {} \n
            time_grad_u23: {} \n
            time_grad_u13: {} \n
            time_grad_u12: {} \n
            time_grad_total: {} \n
        """.format(time_grad_u11_e - time_grad_u11_s,
                   time_grad_u22_e - time_grad_u22_s,
                   time_grad_u33_e - time_grad_u33_s,
                   time_grad_u23_e - time_grad_u23_s,
                   time_grad_u13_e - time_grad_u13_s,
                   time_grad_u12_e - time_grad_u12_s,
                   self.time_total_grad))
        #=======================================================================
        print('Computational time - homogenized coefficient: ')
        self.time_total_assemble = (time_A1212_e - time_A1111_s)
        print("""
            time_A1111: {} \n
            time_A2211: {} \n
            time_A1122: {} \n
            time_A2222: {} \n
            time_A1212: {} \n
            time_total_assemble: {} \n
        """.format(time_A1111_e - time_A1111_s,
                    time_A2222_e - time_A2222_s,
                    time_A3333_e - time_A3333_s,
                    #===========================================================
                    time_A1122_e - time_A1122_s,
                    time_A2211_e - time_A2211_s,
                    time_A1133_e - time_A1133_s,
                    time_A3311_e - time_A3311_s,
                    time_A2233_e - time_A2233_s,
                    time_A3322_e - time_A3322_s,
                    #===========================================================
                    time_A2323_e - time_A2323_s, 
                    time_A1313_e - time_A1313_s, 
                    time_A1212_e - time_A1212_s, 
                    self.time_total_assemble))

        #=======================================================================
        # The homogenized tensor 3-D
        A_ij = np.array([
            [A1111_assembled, A1122_assembled, A1133_assembled, 0, 0, 0],
            [A2211_assembled, A2222_assembled, A2233_assembled, 0, 0, 0],
            [A3311_assembled, A3322_assembled, A3333_assembled, 0, 0, 0],
            [0, 0, 0, A2323_assembled, 0, 0],
            [0, 0, 0, 0, A1313_assembled, 0],
            [0, 0, 0, 0, 0, A1212_assembled],
        ])
        print("The homogenized effective coefficient tensor: \n", np.around(A_ij,1) )

        # ----------------------------------------------------------------------
        """ The homogenized material properties"""
        # [1] Penta, R. and Gerisch, A., 2017. The asymptotic homogenization 
        # elasticity tensor properties for composites with material discontinuities.
        # Continuum Mechanics and Thermodynamics, 29, pp.187-206.

        #=======================================================================
        # [2] Slawinski, M.A., 2010. Waves and rays in elastic continua. World
        #  Scientific. 
        # 5.7. Orthotropic continuum
        # https://web.archive.org/web/20090210192845/http://samizdat.mines.edu/wavesandrays/WavesAndRays.pdf
        # [3] Oliveira, J.A., Pinho-da-Cruz, J. and Teixeira-Dias, F., 2009. 
        # Asymptotic homogenisation in linear elasticity. Part II: Finite 
        # element procedures and multiscale applications. Computational 
        # Materials Science, 45(4), pp.1081-1096.
        #=======================================================================
        
        inv_A_ij = np.linalg.inv(A_ij)
        # print('The inverse homogenized effective coefficient tensor (GPa): \n', inv_A_ij / 1e3)
        E_h_1 = 1/inv_A_ij[0, 0] / 1e3
        E_h_2 = 1/inv_A_ij[1, 1] / 1e3
        E_h_3 = 1/inv_A_ij[2, 2] / 1e3
        G_h_23 = 1/(inv_A_ij[3, 3]) / 1e3
        G_h_13 = 1/(inv_A_ij[4, 4]) / 1e3
        G_h_12 = 1/(inv_A_ij[5, 5]) / 1e3
        nu_h_23 = - inv_A_ij[2, 1] * (1/inv_A_ij[1, 1])
        nu_h_13 = - inv_A_ij[2, 0] * (1/inv_A_ij[0, 0])
        nu_h_12 = - inv_A_ij[1, 0] * (1/inv_A_ij[0, 0])

        print('E_h_1 (GPa): ', round(E_h_1,1))
        print('E_h_2 (GPa): ', round(E_h_2,1))
        print('E_h_3 (GPa): ', round(E_h_3,1))
        # print('G_h_23 (GPa): ', G_h_23)
        # print('G_h_13 (GPa): ', G_h_13)
        # print('G_h_12 (GPa): ', G_h_12)
        print('nu_h_23 (GPa): ', round(nu_h_23,3))
        print('nu_h_13 (GPa): ', round(nu_h_13,3))
        print('nu_h_12 (GPa): ', round(nu_h_12,3))

        material_properties_homogenized = [E_h_1, E_h_2, E_h_3, nu_h_23, nu_h_13, nu_h_12]
        # print(f"E_hom: {E_hom} (MPa)")
        # print(f"nu_hom: {nu_hom}")
        return A_ij, material_properties_homogenized