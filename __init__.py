try:
    from .FEniCS_plot import FEniCS_plot_mode
except:
    print('Failed to import FEniCS_plot')

try:
    from .FEniCS_save_load_solution import \
        Save_XDMF, Save_HDF5, \
        Load_HDF5, Load_XDMF
except:
    print('Failed to import FEniCS_save_load_solution')

# try:
#     from .FEniCS_nonmatching_interpolation import \
#         test_functional2D, test_functional3D
# except:
#     print('Failed to import FEniCS_nonmatching_interpolation')

# try:
#     from .FEniCS_lagrange_interpolator import \
#         test_functional2D, test_functional3D
# except:
#     print('Failed to import FEniCS_lagrange_interpolator')

try:
    from .FEniCS_read_XDMFFile import Read_RB_XDMFFile
except:
    print('Failed to import Read_XDMFFile')
