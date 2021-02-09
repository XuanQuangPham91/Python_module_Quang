# ------------------------------------------------------------------------------
""" FEniCS auxiliary funcions """
# ------------------------------------------------------------------------------

try:
    from .FEniCS_module.FEniCS_plot import (FEniCS_plot_mode)
except:
    print('Failed to import FEniCS_plot')

try:
    from .FEniCS_module.FEniCS_save_load_solution import (
        # save file
        save_XDMF,
        save_HDF5,
        save_vtk,
        save_u_txt,
        # load file
        load_HDF5,
        load_XDMF)
except:
    print('Failed to import FEniCS_save_load_solution')

# try:
#     from .FEniCS_module.FEniCS_lagrange_interpolator import \
#         test_functional2D, test_functional3D
# except:
#     print('Failed to import FEniCS_lagrange_interpolator')

try:
    from .FEniCS_module.FEniCS_read_XDMFFile import (read_RB_XDMFFile)
except:
    print('Failed to import FEniCS_read_XDMFFile')

try:
    from .FEniCS_module.FEniCS_post_processing import (
        # post-processing
        normalize_solution,
        absolute_error,
        relative_error,
        cal_u_magnitude,
        cal_von_mises_stress)
except:
    print('Failed to import FEniCS_auxiliary_code')

# ------------------------------------------------------------------------------
""" Quang auxiliary funcions """
# ------------------------------------------------------------------------------
try:
    from .Quang_module.plot import (
        # plotting
        # Plot_rel_error_norm,
        plot_convergence,
        # plot_general  # has been replace by plot_XY
    )
except:
    print('Failed to import Quang_module.plot')

try:
    from .Quang_module.comparison import (
        # plotting
        plot_XY,
        plot_comparison)
except:
    print('Failed to import Quang_module.comparison')