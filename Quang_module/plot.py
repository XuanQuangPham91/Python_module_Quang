import matplotlib.pyplot as plt
# from dolfin import *
from Python_module_Quang import *
import matplotlib.pyplot as plt
import matplotlib as mpl
# from numpy.core.shape_base import block
from numpy.lib.npyio import loadtxt
# import numpy as np
# import math

sty = 'default'  # "seaborn", "default"
mpl.style.use(sty)
format = 'png'
N = 96

# ------------------------------------------------------------------------------
### Unique functions ###
# ------------------------------------------------------------------------------


def Plot_rel_error_norm(non_affine_rel_errornorm_list,
                        affine_rel_errornorm_list, R_list, fig_path):
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))
    ax1.plot(R_list,
             non_affine_rel_errornorm_list,
             '-ks',
             linewidth=0.8,
             markersize=4.5,
             label='non-affine')
    ax1.plot(R_list,
             affine_rel_errornorm_list,
             '--ko',
             linewidth=0.8,
             markersize=4.5,
             label='affine')
    ax1.set_xlabel('Radius $r$', fontsize=14)
    ax1.set_ylabel(r'Relative error norm $\eta$ [\%]', fontsize=14)
    # ax1.set_yscale('log')
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
    ax1.set_xticks(R_list)
    ax1.grid(color='black', linestyle='--', linewidth=0.5)
    ax1.set_xticks(np.linspace(0.8, 0.9, 11))
    ax1.tick_params(labelsize=12)
    ax1.legend(fontsize=12)
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    plt.close()


def plot_convergence(X, Y, format='png', fig_path=None):
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))
    ax1.plot(X, Y, '-ks', linewidth=0.8, markersize=4.5, label='convergence')
    ax1.set_xlabel('Number of nodes', fontsize=11)
    ax1.set_ylabel('Max of $u_{magnitude}$', fontsize=11)
    # ax1.set_yscale('log')
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
    # ax1.set_xticks(X)
    ax1.grid(color='black', linestyle='--', linewidth=0.5)
    # ax1.set_xticks(Y)
    ax1.tick_params(labelsize=12)
    ax1.legend(fontsize=12)
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    scale_pow = 2
    for x, y in zip(X, Y):
        label = "{:.2f}".format(y * 10**scale_pow)
        plt.annotate(
            label,  # this is the text
            (x, y),  # this is the point to label
            textcoords="offset points",  # how to position the text
            xytext=(0, 5),  # distance from text to points (x,y)
            ha='center')  # horizontal alignment can be left, right or center
    plt.title("Convergence Study", fontsize=13)
    plt.savefig('solution/%s.%s' % ("Convergence_study", format),
                format='%s' % format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    plt.close()


# def plot_convergence(X, Y, fig_path=None, format='png'):
#     fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))
#     ax1.plot(X, Y, '-ks', linewidth=0.8, markersize=4.5, label='convergence')
#     ax1.set_xlabel('Number of nodes', fontsize=11)
#     ax1.set_ylabel('Max of $u_{magnitude}$', fontsize=11)
#     # ax1.set_yscale('log')
#     # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
#     # ax1.set_xticks(X)
#     ax1.grid(color='black', linestyle='--', linewidth=0.5)
#     # ax1.set_xticks(Y)
#     ax1.tick_params(labelsize=12)
#     ax1.legend(fontsize=12)
#     ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
#     ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
#     scale_pow = 2
#     for x, y in zip(X, Y):
#         label = "{:.2f}".format(y * 10**scale_pow)
#         plt.annotate(
#             label,  # this is the text
#             (x, y),  # this is the point to label
#             textcoords="offset points",  # how to position the text
#             xytext=(0, 5),  # distance from text to points (x,y)
#             ha='center')  # horizontal alignment can be left, right or center
#     plt.title("Convergence Study", fontsize=13)
#     plt.savefig('solution/%s.%s' % ("Convergence_study", format),
#                 format='%s' % format)
#     # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
#     plt.close()


def plot_general(X, Y, xlabel, ylabel, title, format='png'):
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    plt.jet()
    ax1.plot(X, Y, '-ks', linewidth=0.8, markersize=4.5, label=title)
    ax1.set_xlabel(xlabel, fontsize=11)
    ax1.set_ylabel(ylabel, fontsize=11)
    # ax1.set_yscale('log')
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
    # ax1.set_xticks(X)
    ax1.grid(color='black', linestyle='--', linewidth=0.5)
    # ax1.set_xticks(Y)
    ax1.tick_params(labelsize=12)
    ax1.legend(fontsize=12)
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    scale_pow = 0
    for x, y in zip(X, Y):
        label = "{:.2f}".format(y * 10**scale_pow)
        plt.annotate(
            label,  # this is the text
            (x, y),  # this is the point to label
            textcoords="offset points",  # how to position the text
            xytext=(0, 5),  # distance from text to points (x,y)
            ha='center')  # horizontal alignment can be left, right or center
    plt.title(title, fontsize=13)
    plt.savefig('solution/%s.%s' % (title, format), format='%s' % format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    plt.close()


def comparison_various_mesh(X=None, Y=None):
    path32 = "solution/comparison_80/Relative_errornorm_32.txt"
    path80 = "solution/comparison_80/Relative_errornorm_80.txt"
    path96 = "solution/comparison_80/Relative_errornorm_96.txt"
    Y = loadtxt(path32, delimiter=',')
    Y2 = loadtxt(path80, delimiter=',')
    Y96 = loadtxt(path96, delimiter=',')

    X = [
        '$\mu_0$ \n (0.5, 0.5)', '$\mu_1$ \n (0.75, 0.75)',
        '$\mu_2$\n (1.0, 1.0)', '$\mu_3$\n (1.25, 1.25)',
        '$\mu_4$\n (1.5, 1.5)'
    ]

    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    format = 'png'
    title = ("Relative errornorm")
    ax1.plot(X, Y, '-ks', linewidth=0.8, markersize=4.5, label=title + " 32")
    ax1.plot(X, Y2, '-bs', linewidth=0.8, markersize=4.5, label=title + " 80")
    ax1.plot(X, Y96, '-gs', linewidth=0.8, markersize=4.5, label=title + " 96")
    xlabel = "$\mu$"
    ylabel = "Relative errornorm ($\%$)"
    ax1.set_xlabel(xlabel, fontsize=11)
    ax1.set_ylabel(ylabel, fontsize=11)
    # ax1.set_yscale('log')
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
    # ax1.set_xticks(X)
    ax1.grid(color='black', linestyle='--', linewidth=0.5)
    # ax1.set_xticks(Y)
    ax1.tick_params(labelsize=12)
    ax1.legend(fontsize=12)
    # ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    scale_pow = 0
    # for x, y in zip(X, Y):
    #     label = "{:.2f}".format(y * 10**scale_pow)
    #     plt.annotate(
    #         label,  # this is the text
    #         (x, y),  # this is the point to label
    #         textcoords="offset points",  # how to position the text
    #         xytext=(0, 5),  # distance from text to points (x,y)
    #         ha='center')  # horizontal alignment can be left, right or center
    # plt.title(title, fontsize=13)
    # plt.savefig('solution/%s.%s' % (title, format), format='%s' % format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    # plt.close()


def comparison_various_80():
    # variable list = Y80_0, Y80_1, Y80_2
    # for index in range(4):
    #     # exec(f'u_{index} = k')or index in range(5):
    #     # exec(f'u_{index} = k')
    #     Y = loadtxt("solution/comparison_80/Relative_errornorm_80_%d.txt" %
    #                 index,
    #                 delimiter=',')
    #     exec(f'Y80_{index} = Y')
    # title = 'Y80_%s' % index
    path0 = "solution/comparison_80/Relative_errornorm_80_0.txt"
    path1 = "solution/comparison_80/Relative_errornorm_80_0.txt"
    path2 = "solution/comparison_80/Relative_errornorm_80_0.txt"
    path3 = "solution/comparison_80/Relative_errornorm_80_0.txt"
    Y80_0 = loadtxt(path0, delimiter=',')
    Y80_1 = loadtxt(path1, delimiter=',')
    Y80_2 = loadtxt(path2, delimiter=',')
    Y80_3 = loadtxt(path3, delimiter=',')

    X = [
        '$\mu_0$ \n (0.5, 0.5)', '$\mu_1$ \n (0.75, 0.75)',
        '$\mu_2$\n (1.0, 1.0)', '$\mu_3$\n (1.25, 1.25)',
        '$\mu_4$\n (1.5, 1.5)'
    ]

    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    format = 'png'
    title = ("Relative errornorm")
    ax1.set_title(
        ("Re-calculate relative errornorm various times \n with same number of nodes"
         ).format(sty),
        fontsize=13,
        # color='C0'
    )
    ax1.plot(X,
             Y80_0,
             '-rs',
             linewidth=0.8,
             markersize=4.5,
             label=title + " 0")
    ax1.plot(X,
             Y80_1,
             '-bs',
             linewidth=0.8,
             markersize=4.5,
             label=title + " 1")
    ax1.plot(X,
             Y80_2,
             '-gs',
             linewidth=0.8,
             markersize=4.5,
             label=title + " 2")
    ax1.plot(X,
             Y80_3,
             '-ks',
             linewidth=0.8,
             markersize=4.5,
             label=title + " 3")
    xlabel = "$\mu$"
    ylabel = "Relative errornorm ($\%$)"
    ax1.set_xlabel(xlabel, fontsize=11)
    ax1.set_ylabel(ylabel, fontsize=11)
    # ax1.set_yscale('log')
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
    # ax1.set_xticks(X)
    ax1.grid(color='black', linestyle='--', linewidth=0.5)
    # ax1.set_Yticks(Y)
    ax1.tick_params(labelsize=12)
    ax1.legend(fontsize=11)
    # ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # scale_pow = 0
    # for x, y in zip(X, Y):
    #     label = "{:.2f}".format(y * 10**scale_pow)
    #     plt.annotate(
    #         label,  # this is the text
    #         (x, y),  # this is the point to label
    #         textcoords="offset points",  # how to position the text
    #         xytext=(0, 5),  # distance from text to points (x,y)
    #         ha='center')  # horizontal alignment can be left, right or center
    # plt.title(title, fontsize=13)
    plt.savefig('solution/%s.%s' % ((title + " various times"), format),
                format='%s' % format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    # plt.close()


# ------------------------------------------------------------------------------
### General functions ###
# ------------------------------------------------------------------------------


def plot_comparison(title,
                    titleY1,
                    titleY2,
                    format='png',
                    xlabel="x",
                    ylabel="y",
                    X=None,
                    Y1=None,
                    Y2=None,
                    pathX=None,
                    pathY1=None,
                    pathY2=None,
                    style=None):
    # path1 = "solution/comparison_80/Relative_errornorm_32.txt"
    # path2 = "solution/comparison_80/Relative_errornorm_80.txt"

    if pathX != None:
        X = loadtxt(pathX, delimiter=',')
    if pathY1 != None:
        Y1 = loadtxt(pathY1, delimiter=',')
    if pathY1 != None:
        Y2 = loadtxt(pathY2, delimiter=',')
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    title = ("Relative errornorm")
    ax1.set_title(
        title.format(sty),
        fontsize=13,
        # color='C0'
    )
    ax1.plot(X, Y1, '-ks', linewidth=0.8, markersize=4.5, label=titleY1)
    ax1.plot(X, Y2, '-bs', linewidth=0.8, markersize=4.5, label=titleY2)
    ax1.set_xlabel(xlabel, fontsize=11)
    ax1.set_ylabel(ylabel, fontsize=11)
    # ax1.set_yscale('log')
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
    # ax1.set_xticks(X)
    ax1.grid(color='black', linestyle='--', linewidth=0.5)
    # ax1.set_xticks(Y)
    ax1.tick_params(labelsize=11)
    ax1.legend(fontsize=11)
    if style == 'sci':
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # scale_pow = 0

    ### Data label ###
    # for x, y in zip(X, Y):
    #     label = "{:.2f}".format(y * 10**scale_pow)
    #     plt.annotate(
    #         label,  # this is the text
    #         (x, y),  # this is the point to label
    #         textcoords="offset points",  # how to position the text
    #         xytext=(0, 5),  # distance from text to points (x,y)
    #         ha='center')  # horizontal alignment can be left, right or center

    ### save figure ###
    plt.savefig('solution/%s.%s' % (title, format), format=format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    plt.close()


def plot_XY(title,
            format='png',
            xlabel="x",
            ylabel="y",
            X=None,
            Y=None,
            pathX=None,
            pathY=None,
            style=None):
    # path1 = "solution/comparison_80/Relative_errornorm_32.txt"  # or csv
    if pathX != None:
        X = loadtxt(pathX, delimiter=',')
    if pathY != None:
        Y = loadtxt(pathY, delimiter=',')
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    ax1.set_title(
        title.format(sty),
        fontsize=13,
        # color='C0'
    )
    ax1.plot(X, Y, '-ks', linewidth=0.8, markersize=4.5, label=title)
    ax1.set_xlabel(xlabel, fontsize=11)
    ax1.set_ylabel(ylabel, fontsize=11)
    # ax1.set_yscale('log')
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 5e0])
    # ax1.set_xticks(X)
    ax1.grid(color='black', linestyle='--', linewidth=0.5)
    # ax1.set_xticks(Y)
    ax1.tick_params(labelsize=11)
    # ax1.legend(fontsize=12)
    if style == 'sci':
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    ## Data label ###
    scale_pow = 0
    for x, y in zip(X, Y):
        label = "{:.2f}".format(y * 10**scale_pow)
        plt.annotate(
            label,  # this is the text
            (x, y),  # this is the point to label
            textcoords="offset points",  # how to position the text
            xytext=(0, 5),  # distance from text to points (x,y)
            ha='center')  # horizontal alignment can be left, right or center

    ### save figure ###
    plt.savefig('solution/%s.%s' % (title, format), format=format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    # plt.close()


# ------------------------------------------------------------------------------
### Main program ###
# ------------------------------------------------------------------------------

# if __name__ == "__main__":
#     comparison_various_mesh()
#     comparison_various_80()

#     plt.show(block=False)
