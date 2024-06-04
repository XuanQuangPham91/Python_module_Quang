# ==============================================================================
""" FEniCS auxiliary funcions """
# ==============================================================================
try:
    from Python_module_Quang.FEniCS_module.FEniCS_plot import (FEniCS_plot_mode)
except:
    print('Failed to import FEniCS_plot')

try:
    from Python_module_Quang.FEniCS_module.FEniCS_save_load_solution import (
        ### save file ##
        save_XDMF,
        save_HDF5,
        save_vtk,
        save_u_txt,
        ### load file ###
        load_HDF5,
        load_XDMF,
        load_mesh_RBniCS,
        RBniCS_convert_mesh)
except:
    print('Failed to import FEniCS_save_load_solution')

# try:
#     from .FEniCS_module.FEniCS_lagrange_interpolator import \
#         test_functional2D, test_functional3D
# except:
#     print('Failed to import FEniCS_lagrange_interpolator')

try:
    from Python_module_Quang.FEniCS_module.FEniCS_read_XDMFFile import (read_RB_XDMFFile)
except:
    print('Failed to import FEniCS_read_XDMFFile')

try:
    from Python_module_Quang.FEniCS_module.FEniCS_post_processing import (
        ### post-processing ###
        normalize_solution,
        absolute_error,
        relative_error,
        cal_u_magnitude,
        cal_von_mises_stress)
except:
    print('Failed to import FEniCS_auxiliary_code')

# ==============================================================================
""" Quang auxiliary funcions """
# ==============================================================================
try:
    from Python_module_Quang.Quang_module.plot import (
        ### plotting ###
        # Plot_rel_error_norm,
        plot_convergence,
        plot_XY,
        plot_comparison,
        # plot_general  # has been replace by plot_XY
    )
except:
    print('Failed to import Quang_module.plot')

# ==============================================================================
""" Homogenization_module """
# ==============================================================================

try:
    from Python_module_Quang.Homogenization_module.auxiliary_codes import (
        # Function =============================================================
        check_create_dir,
        delete_model,
        check_log,
        copy_and_overwrite,
        RBniCS_convert_mesh,
        my_plot,
        plot_greedy,
        plot_greedy_deim,
        plot_sample_with_error_in_effecitivity,
        plot_box_whisker,
        plot_error_analysis,
        plot_error_analysis_manual,
        plot_error_analysis_effectivity,
        manual_error_analysis,
        basis_function_grad,
        relocation_offset_text,
        # Class ================================================================
        Logger,
        PeriodicBoundary,
        PeriodicBoundary3D,
        error_estimation,
        Pinpoint,
    )
except:
    print('Failed to import Homogenization_module.auxiliary_codes')


try:
    from Python_module_Quang.Homogenization_module.LatinHyperCube_sampling import (
        # Function =============================================================
        sample_set_RBniCS_formater,
        _Latin_Hypercube_sampling,
        Latin_Hypercube_sampling_func,
        # Class ================================================================
    )
except:
    print('Failed to import Homogenization_module.LatinHyperCube_sampling')

try:
    from Python_module_Quang.Homogenization_module.homogenizedTensor import (
        # Function =============================================================
        # Class ================================================================
        homogenizedTensor,
        homogenizedTensor3D
    )
except:
    print('Failed to import Homogenization_module.homogenizedTensor')

try:
    from Python_module_Quang.Gmsh_container.mesh_converter_ASCII_v2_PRIOR import (
        # Function =============================================================
        gmsh_xml_converter_2D
        # Class ================================================================
    )
except:
    print('Failed to import Homogenization_module.homogenizedTensor')