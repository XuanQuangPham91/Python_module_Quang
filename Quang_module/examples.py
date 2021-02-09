import matplotlib.pyplot as plt


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


def plot_convergence(X, Y, fig_path=None):
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
    format = 'png'
    plt.savefig('solution/%s.%s' % ("Convergence_study", format),
                format='%s' % format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    plt.close()


def plot_convergence(X, Y, fig_path=None):
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
    format = 'png'
    plt.savefig('solution/%s.%s' % ("Convergence_study", format),
                format='%s' % format)
    # fig.savefig(fig_path, bbox_inches='tight', dpi=400)
    plt.close()


def plot_general(X, Y, xlabel, ylabel, title, format='png'):
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
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
    # plt.close()
