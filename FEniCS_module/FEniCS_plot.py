# Define auxiliary code
''' visualization  '''
from dolfin import plot
import matplotlib.pyplot as plt
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


def FEniCS_plot_mode(u,
                     title,
                     mode="glyphs",
                     number_of_figure=None,
                     reduced_problem=None,
                     savefig=None,
                     grid=None):
    if number_of_figure == None:
        pass
    elif number_of_figure == True:
        plt.figure()
    else:
        plt.figure(number_of_figure)
    if reduced_problem == 'reduced_problem':
        # u = reduced_problem
        p = plot(u, reduced_problem=reduced_problem, title=title, mode=mode)
    elif reduced_problem == None:
        p = plot(u, title='%s' % title, mode=mode)
    # plt.axis('tight')
    # plt.legend()
    if grid == True:
        plt.grid(True)
    else:
        pass
    plt.title(title)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.colorbar(p)
    if savefig == True:
        plt.savefig('solution/%s.png' % title, format='png')
    else:
        pass


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