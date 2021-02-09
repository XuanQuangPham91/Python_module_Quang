# from dolfin import *
# from rbnics import *
from Python_module_Quang import *
# import clear_offline_data
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy.core.shape_base import block
from numpy.lib.npyio import loadtxt
# import numpy as np
# import math

sty = 'default'  # "seaborn", "default"
mpl.style.use(sty)
format = 'png'
N = 96


def comparison_various_mesh():
    Y32 = loadtxt("solution/comparison_32/Relative_errornorm_32.txt",
                  delimiter=',')
    Y80 = loadtxt("solution/comparison_80/Relative_errornorm_80_0.txt",
                  delimiter=',')
    Y96 = loadtxt("solution/comparison_96/Relative_errornorm_96.txt",
                  delimiter=',')

    print(Y32)
    # exit()

    X = [
        '$\mu_0$ \n (0.5, 0.5)', '$\mu_1$ \n (0.75, 0.75)',
        '$\mu_2$\n (1.0, 1.0)', '$\mu_3$\n (1.25, 1.25)',
        '$\mu_4$\n (1.5, 1.5)'
    ]

    xlabel = "$\mu$"
    ylabel = "Relative errornorm ($\%$)"
    title = ("Relative errornorm")
    format = 'png'
    # plt.figure()

    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    # plt.jet()
    ax1.plot(X, Y32, '-ks', linewidth=0.8, markersize=4.5, label=title + " 32")
    ax1.plot(X, Y80, '-bs', linewidth=0.8, markersize=4.5, label=title + " 80")
    ax1.plot(X, Y96, '-gs', linewidth=0.8, markersize=4.5, label=title + " 96")
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
    # for x, y in zip(X, Y32):
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

    Y80_0 = loadtxt("solution/comparison_80/Relative_errornorm_80_0.txt",
                    delimiter=',')
    Y80_1 = loadtxt("solution/comparison_80/Relative_errornorm_80_1.txt",
                    delimiter=',')
    Y80_2 = loadtxt("solution/comparison_80/Relative_errornorm_80_2.txt",
                    delimiter=',')
    Y80_3 = loadtxt("solution/comparison_80/Relative_errornorm_80_3.txt",
                    delimiter=',')

    X = [
        '$\mu_0$ \n (0.5, 0.5)', '$\mu_1$ \n (0.75, 0.75)',
        '$\mu_2$\n (1.0, 1.0)', '$\mu_3$\n (1.25, 1.25)',
        '$\mu_4$\n (1.5, 1.5)'
    ]

    xlabel = "$\mu$"
    ylabel = "Relative errornorm ($\%$)"
    title = ("Relative errornorm")
    format = 'png'
    # plt.figure()

    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    # plt.jet()
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
    # for x, y in zip(X, Y32):
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


if __name__ == "__main__":
    # comparison_various_mesh()
    comparison_various_80()

    plt.show(block=False)
