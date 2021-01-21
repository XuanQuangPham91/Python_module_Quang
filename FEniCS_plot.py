# Define auxiliary code
''' visualization  '''
import matplotlib.pyplot as plt
import os

print(os.getcwd())


def my_plot(
    number_of_figure,
    u,
    title,
    reduced_problem='False',
    savefig="True",
):
    if reduced_problem == 'reduced_problem':
        # u = reduced_problem
        p = plot(u, reduced_problem=reduced_problem, title=title)
    elif reduced_problem == 'False':
        p = plot(u, title='%s' % title)
    plt.figure(number_of_figure)
    plt.colorbar(p)
    # plt.axis('tight')
    # plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    if savefig == 'True':
        plt.savefig('20210118_RB_results/%s.png' % title, format='png')
    else:
        pass


def my_plot_mode(number_of_figure,
                 u,
                 title,
                 mode,
                 reduced_problem='False',
                 savefig="True"):
    if reduced_problem == 'reduced_problem':
        # u = reduced_problem
        p = plot(u, reduced_problem=reduced_problem, title=title, mode=mode)
    elif reduced_problem == 'False':
        p = plot(u, title='%s' % title, mode=mode)
    plt.figure(number_of_figure)
    # plt.axis('tight')
    # plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.colorbar(p)
    if savefig == 'True':
        plt.savefig('20210118_RB_results/%s.png' % title, format='png')
    else:
        pass