# Define auxiliary code
''' visualization  '''
from dolfin import plot
import matplotlib.pyplot as plt
from numpy.lib.shape_base import tile
from .FEniCS_save_load_solution import save_txt
# import os

# print(os.getcwd())

#------------------------------------------------------------------------------
'''FEniCS plot '''

# def FEniCS_plot_noneMode(
#     number_of_figure,
#     u,
#     title,
#     reduced_problem=None,
#     savefig=None,
# ):
#     if reduced_problem == 'reduced_problem':
#         # u = reduced_problem
#         p = plot(u, reduced_problem=reduced_problem, title=title)
#     elif reduced_problem == None:
#         p = plot(u, title='%s' % title)
#     plt.figure(number_of_figure)
#     plt.colorbar(p)
#     # plt.axis('tight')
#     # plt.legend()
#     plt.grid(True)
#     plt.title(title)
#     plt.xlabel('$x$')
#     plt.ylabel('$y$')
#     if savefig == True:
#         plt.savefig('solution/%s.png' % title, format='png')
#     else:
#         pass


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
        p = plot(u, title='%s' % title)
    else:
        p = plot(u, title='%s' % title, mode=mode)
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
        plt.savefig('solution/%s.%s' % (title, format), format='%s' % format)
    else:
        pass
    if savefiletxt == True:
        save_txt(u, title=title)
    else:
        pass
    plt.close()


# def my_FEniCS_plot(number_of_figure, u, title):
#     plt.figure(number_of_figure)
#     plt.colorbar(plot(u, title='%s' % title))
#     # plt.axis('tight')
#     # plt.legend()
#     # plt.grid(True)
#     plt.title(title)
#     plt.xlabel('$y_1$')
#     plt.ylabel('$y_2$')
#     plt.savefig('solution/%s.png' % title, format='png')

# def my_FEniCS_plot_mode(number_of_figure, u, title, mode):
#     plt.figure(number_of_figure)
#     plt.colorbar(plot(u, title='%s' % title, mode=mode))
#     # plt.axis('tight')
#     # plt.legend()
#     # plt.grid(True)
#     plt.title(title)
#     plt.xlabel('$y_1$')
#     plt.ylabel('$y_2$')
#     plt.savefig('solution/%s.png' % title, format='png')