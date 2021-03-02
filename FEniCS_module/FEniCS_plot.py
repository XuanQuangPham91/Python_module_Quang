# Define auxiliary code
''' visualization  '''
from dolfin import plot
import matplotlib.pyplot as plt
# from numpy.lib.shape_base import tile
from .FEniCS_save_load_solution import save_u_txt
# import os

# print(os.getcwd())

#------------------------------------------------------------------------------
'''FEniCS plot '''
#------------------------------------------------------------------------------


def FEniCS_plot_mode(
    u,
    title,
    format="png",
    xlabel='$x$',
    ylabel='$y$',
    mode=None,
    number_of_figure=True,
    savefig=None,
    savefiletxt=None,
    grid=None,
):
    """
    mode=None, "glyphs", "displacement",
    """
    plt.jet()
    if number_of_figure == None:
        pass
    elif number_of_figure == True:
        plt.figure()
    else:
        plt.figure(number_of_figure)
    if mode == None:
        p = plot(u, title=title)
    else:
        p = plot(u, title=title, mode=mode)
    # plt.axis('tight')
    # plt.legend()
    if grid == True:
        plt.grid(True)
    else:
        pass
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(p)
    if savefig == True:
        # plt.savefig('solution/%s.%s' % (title, format), format='%s' % format)
        plt.savefig(f'solution/{title}.{format}', format=format)
    else:
        pass
    if savefiletxt == True:
        save_u_txt(u, title=title)
    else:
        pass
    plt.close()
