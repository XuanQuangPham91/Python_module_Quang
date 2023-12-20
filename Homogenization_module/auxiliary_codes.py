from dolfin import *
import dolfin as df
from mshr import *
from Python_module_Quang import *
from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoMinorLocator
# plticker
# ==============================================================================
# ==============================================================================

import copy
import json
import logging
import meshio
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import pandas as pd
import shutil
import sys
import time

# ==============================================================================
# ==============================================================================
from matplotlib import rc

rc('font', **{'family': 'serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
# matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
# matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsfonts}'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsfonts,amsmath,amssymb}'

fontsize = 21
matplotlib.rc('xtick', labelsize=fontsize)
matplotlib.rc('ytick', labelsize=fontsize)

# Elsevier standard: 
## SmallDouble column: 90mm - 3.53in
## halfpage-single column: 140mm - 5.1in
## fullpage-single column: 190mm - 7.48in
Elsevier_standard = "FullSingleColumn"
if Elsevier_standard == "SmallDoubleColumn":
    # value_height =
    value_width = 3.53
if Elsevier_standard == "HalfSingleColumn":
    value_width = 5.1
if Elsevier_standard == "FullSingleColumn":
    value_width = 7.48

start = time.time()

# ==============================================================================
# ==============================================================================
# In[]


def check_create_dir(data_RB_path):
    if os.path.exists(data_RB_path) == False:
        os.mkdir(data_RB_path)

def delete_model(delete_state=False, modelName=""):
    if delete_state and os.path.isdir(modelName):
        shutil.rmtree(modelName)
        print("Removed folder")

def check_log(log_path, remove=False):
    """
    remove == False (default)
    """
    if remove:
        if os.path.isfile(log_path):
            os.remove(log_path)
    sys.stdout = Logger(log_path)


def copy_and_overwrite(from_path, to_path):
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree(from_path, to_path)


class Logger(object):

    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def __getattnu_f__(self, attr):
        return getattr(self.terminal, attr)

    def write(self, message):
        self.log.write(message)

    def flush(self):
        self.log.flush()
#In[] 
class PeriodicBoundary(SubDomain):
    """
    # Periodic boundary condition definition
        * Origin at center 
    """
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of 
        # the two corners (0, 1) and (1, 0)
        self.L = 1.0
        Move = [0.0, 0.0]  # at center of coordinate system
        self.x1 = [-self.L/2 + Move[0], -self.L/2 + Move[1]]
        self.x2 = [-self.L/2 + Move[0], +self.L/2 + Move[1]]
        self.x3 = [+self.L/2 + Move[0], +self.L/2 + Move[1]]
        self.x4 = [+self.L/2 + Move[0], -self.L/2 + Move[1]]
        return bool((near(x[0], self.x1[0]) or near(x[1], self.x1[1]))
                    and (not ((near(x[0], self.x2[0]) and near(x[1], self.x2[1])) or
                              (near(x[0], self.x4[0]) and near(x[1], self.x4[1]))))
                    and on_boundary)

    def map(self, x, y):
        if near(x[0], self.x3[0]) and near(x[1], self.x3[0]):
            y[0] = x[0] - self.L
            y[1] = x[1] - self.L
        elif near(x[0], self.x3[0]):
            y[0] = x[0] - self.L
            y[1] = x[1]
        elif near(x[1], self.x3[1]):
            # else:  # near(x[1], 1)
            y[0] = x[0]
            y[1] = x[1] - self.L


# In[]


def RBniCS_convert_mesh(solution_path, mesh_path):
    """
    https://pythonlang.dev/repo/nschloe-meshio/
    """
    with meshio.xdmf.TimeSeriesReader(solution_path) as reader:
        points, cells = reader.read_points_cells()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)

    mesh = meshio.Mesh(points=points,
                       cells=cells,
                       point_data=point_data,
                       cell_data=cell_data)
    mesh.write(mesh_path)

    mesh = Mesh(mesh_path)
    File(mesh_path) << mesh

    return mesh


# In[]


def my_plot(
    u,
    path='solution/',
    title='missed title',
    xlabel=r'$x_1$',
    ylabel=r'$x_2$',
    format='pdf',
    mode=None,
):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_figwidth(value_width)
    plt.jet()
    if mode == None:
        fig.colorbar(plot(u))
    else:
        fig.colorbar(df.plot(u, mode=mode))
    # plt.axis('tight')
    # plt.legend()
    ax.grid(True)
    # ax.set_title(title, fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_box_aspect(1)

    fig.savefig(f'{path}.{format}',
                format=format,
                bbox_inches='tight',
                dpi=600)
    plt.show(block=False)
    # time.sleep(2)
    plt.close('all')


def plot_greedy(
    RB_folder,
    fig_path,
    legend=None,
):
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        # figsize=(5, 5),
    )
    fig.set_figwidth(value_width)
    max_N_list_affine_list = []
    for i_folder in range(len(RB_folder)):

        greedy_abs_erronu_f_file_path = os.path.join(
            RB_folder[i_folder], 'post_processing', 'error_estimator_max.txt')
        with open(greedy_abs_erronu_f_file_path, 'r') as f:
            greedy_abs_error = json.loads(f.read())
        greedy_rel_error = []

        for i, error in enumerate(greedy_abs_error):
            if i == 0:
                greedy_rel_error.append(1)
            else:
                greedy_rel_error.append(greedy_abs_error[i] /
                                        greedy_abs_error[0])

        N_list_affine = [i for i in range(len(greedy_rel_error))]
        max_N_list_affine_list.append(max(N_list_affine))
        color_list = ['-ro', '--b^', '-.gs']
        ax.plot(
            N_list_affine,
            greedy_rel_error,
            color_list[i_folder],
            linewidth=1.0,
            markersize=6.5,
            alpha=0.7,
        )
    if legend == None:
        pass
    else:
        ax.legend(legend, ncol=1, bbox_to_anchor=(1.02, 1.02), 
               fontsize=fontsize, fancybox=True, framealpha=0.9)
    ax.set_yscale('log')
    ax.grid(True)
    # ax.grid(color='black', linestyle='--', linewidth=0.5)
    # ax.set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    # ax.set_title(r"{title}".format(title=title), fontsize=fontsize)

    if (N_list_affine[-1] % 2) == 0:
        ax.set_xlim(-0.25, max(max_N_list_affine_list) + 1.25)
    else:
        ax.set_xlim(-0.25, max(max_N_list_affine_list) + 0.25)
    # ax.set_xlim(-0.25, 14.25)

    ax.set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_xlabel(r'number of basis functions $\textit{N}$', fontsize=fontsize)
    ax.set_ylabel(
        r'Relative RB error estimate $\Delta_{\textrm{rel}}^{N}(\mu)$',
        fontsize=fontsize)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax.set_box_aspect(1)
    ax.tick_params(labelsize=fontsize)
    fig.savefig(fig_path, bbox_inches='tight', dpi=600)
    plt.show(block=False)
    plt.close('all')


# return greedy_rel_error

# In[]
def plot_sample_with_error_in_effecitivity(current_directory, RB_folder,
                                           fig_path):
    error_folder = os.path.join(current_directory, RB_folder,
                                'manual_error_analysis')
    error_data_list = os.listdir(error_folder)
    N = len(error_data_list)
    error_file_path = os.path.join(error_folder, '%s.csv' % N)
    error_data = pd.read_csv(error_file_path)
    effectivity_lower_one = error_data[error_data['effectivity'] < 1]
    effectivity_larger_one = error_data[error_data['effectivity'] >= 1]

    fig, ax = plt.subplots(1,1)
    fig.set_figwidth(value_width)
    ax.scatter(effectivity_larger_one['radius'],
               effectivity_larger_one['K_fiber'],
               marker='x',
               c="blue",
               s=100,
               label='Effectivity $\geqslant 1$')
    ax.scatter(effectivity_lower_one['radius'],
               effectivity_lower_one['K_fiber'],
               marker='o',
               c="none",
               edgecolors='red',
               s=100,
               linewidths=2.0,
               label='Effectivity $< 1$')
    font_size = 14
    tick_size = 12
    ax.grid(color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel('radius $r$', fontsize=font_size)
    ax.set_ylabel('$k_{\mathrm{f}}$', fontsize=font_size)
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=tick_size)
    ax.set_aspect('equal')
    fig.savefig(fig_path, bbox_inches='tight', dpi=600)
    plt.close()


def plot_box_whisker(
    current_directory,
    RB_folder,
    fig_path,
    boxplot_data_path,
    configuration=None,
    plot_header='effectivity',
):

    def get_box_plot_data(N_list, bp):
        rows_list = []
        for i in range(len(N_list)):
            dict1 = {}
            dict1['N'] = N_list[i]
            dict1['lower_whisker'] = bp['whiskers'][i * 2].get_ydata()[1]
            dict1['lower_quartile'] = bp['boxes'][i].get_ydata()[1]
            dict1['median'] = bp['medians'][i].get_ydata()[1]
            dict1['upper_quartile'] = bp['boxes'][i].get_ydata()[2]
            dict1['upper_whisker'] = bp['whiskers'][(i * 2) + 1].get_ydata()[1]
            rows_list.append(dict1)
        return pd.DataFrame(rows_list)

    plt.rcParams['text.usetex'] = True

    error_folder = os.path.join(current_directory, RB_folder,
                                'manual_error_analysis')
    error_data_list = os.listdir(error_folder)
    N_list = [i + 1 for i in range(len(error_data_list))]
    dict = {}
    for N in N_list:
        error_file_path = os.path.join(error_folder, '%s.csv' % N)
        error_data = pd.read_csv(error_file_path)
        eff = error_data[plot_header].values.tolist()
        dict['%s' % N] = eff
    fig, ax = plt.subplots(nrows=1,ncols=1)
    fig.set_figwidth(value_width)
    ax.set_box_aspect(1)
    boxplot = ax.boxplot(dict.values(), showfliers=True)
    font_size = fontsize
    tick_size = fontsize
    # ax.set_xticks(N)
    # ax.set_yscale('log')
    ax.grid(color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel(r'Number of basis functions $\textit{N}$',
                  fontsize=font_size)

    ax.tick_params(labelsize=tick_size)
    xleft, yright = ax.get_xlim()
    interval = 3
    ## =====================================================================
    if interval == 2:
        ## interval = 2
        if yright % 2 == 0:
            yright = yright + 2
            # ax.set_xlim(0.0, yright)
        else:
            yright = yright + 1
            # ax.set_xlim(0.0, yright)
        ax.set_xticks(np.arange(0.0, yright, 2.0))
        ax.set_xticklabels(np.arange(0.0, yright, 2.0, dtype=int))
    ## =====================================================================
    elif interval == 3:
        ## interval = 3
        if yright % 3 == 0:
            yright = yright + 1
            # ax.set_xlim(0.0, yright)
        else:
            yright = yright + 1
            while yright % 3 == 0:
                yright = yright + 1
                # ax.set_xlim(0.0, yright)
            # ax.set_xlim(0.0, yright)
        ax.set_xticks(np.arange(0.0, yright, 3.0))
        ax.set_xticklabels(np.arange(0.0, yright, 3.0, dtype=int))
    ## =========================================================================
    ## =========================================================================
    # if yright % 2 == 0:
    #     yright = yright + 2
    #     # ax.set_xlim(0.0, yright)
    # else:
    #     yright = yright + 1
    #     # ax.set_xlim(0.0, yright)
    # ax.set_xticks(np.arange(0.0, yright, 2.0))
    # ax.set_xticklabels(np.arange(0.0, yright, 2.0, dtype=int))
    ## =========================================================================
    ## =========================================================================

    # locator = MultipleLocator(10)
    # ax.xaxis.set_major_locator(locator)
    # ax.xaxis.set_minor_locator(AutoMinorLocator())
    # loc = plticker.MultipleLocator(
    #     base=2.0)  # this locator puts ticks at regular intervals
    # ax.xaxis.set_major_locator(loc)

    if plot_header == 'effectivity':
        ax.set_ylabel(r'Effectivity $\eta^{N}_\mathbb{V}(\mu)$', fontsize=font_size)
        # ax.set_yticks(range(int(np.round(ymax, 0)) + 1))
        if configuration == "Oliveira":
            ax.set_ylim(0, 16)  # oliveira
        if configuration == "Macedo":
            ax.set_ylim(0, 45)  # macedo
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    elif plot_header == 'residual_norm':
        # ax.set_ylabel(r'Residual norm $\lVert \hat{r}_\delta (\mu) \rVert$',
        ax.set_ylabel(r'Dual norm of residual $ \| \hat{r}_\delta (\mu) \|$',
                      fontsize=font_size)
        ax.set_yscale('log')
    elif plot_header == 'residual_norm_squared':
        ax.set_ylabel(
            r'Squared of dual norm of residual $\| \hat{r}_\delta (\mu) \| ^2$',
            fontsize=font_size)
        ax.set_yscale('log')
    # elif plot_header == 'SCM_alpha_LB':
    #     ax.set_ylabel(r'Truth error $ \| e (\mu) \|_{\mathbb{V}} $',
    #                   fontsize=font_size)
    elif plot_header == 'error':  #\mathrm{\mathbb{}} ## error exact
        ax.set_ylabel(r'Truth error $ \| e (\mu) \|_{\mathbb{V}} $',
                      fontsize=font_size)
        ax.set_yscale('log')
        # ax.set_yticklabels([])
    elif plot_header == 'error_estimator':
        ax.set_ylabel(r'Error estimator $\Delta_{N} (\mu)$',
                      fontsize=font_size)
        ax.set_yscale('log')
        # ax.set_yticklabels([])

    # ax.set_ylim((0,12))
    # ax.set_yticks([0, 1, 2, 3, 4, 5])
    # ax.set_yticks(range(int(np.round(ymax, 0)) + 1))
    # ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_yticks([0,5,10,15,20])
    # ax.legend(fontsize=font_size)
    fig.savefig(fig_path, bbox_inches='tight', dpi=600)
    boxplot_data = get_box_plot_data(N_list, boxplot)
    boxplot_data.to_csv(boxplot_data_path, index=False)
    plt.close()


# In[]


def plot_error_analysis(current_directory, RB_folder, fig_path):
    error_file_path = os.path.join(current_directory, RB_folder,
                                   'error_analysis', 'error_analysis',
                                   'solution_u_error.csv')
    error_data = pd.read_csv(error_file_path, delimiter=';')
    N = error_data['N'].values.tolist()
    # print(f"N: {N}")
    max_error = error_data['max(error_u)'].values.tolist()
    # print(f"max_error: {max_error}")
    max_bound = error_data['max(error_estimator_u)'].values.tolist()
    # print(f"max_bound: {max_bound}")
    mean_error = error_data['gmean(error_u)'].values.tolist()
    # print(f"mean_error: {mean_error}")
    mean_bound = error_data['gmean(error_estimator_u)'].values.tolist()
    fig, ax = plt.subplots(nrows=1,ncols=1)
    fig.set_figwidth(value_width)
    ax.plot(N,
            max_bound,
            ':bo',
            linewidth=3.0,
            markersize=6.5,
            label=r'Max $\Delta^{N}$')
    ax.plot(N,
            max_error,
            '-bo',
            linewidth=1.0,
            markersize=8.5,
            label=r'Max $||e||_X$',
            mfc='none')
    ax.plot(N,
            mean_bound,
            ':r^',
            linewidth=3.0,
            markersize=6.5,
            label=r'Average $\Delta^{N}$')
    ax.plot(N,
            mean_error,
            '-r^',
            linewidth=1.0,
            markersize=8.5,
            label=r'Average $||e||_X$',
            mfc='none')

    font_size = fontsize
    tick_size = fontsize
    # ax.set_xticks(N)
    xleft, yright = ax.get_xlim()
    interval = 3
    ## =====================================================================
    if interval == 2:
        ## interval = 2
        if yright % 2 == 0:
            yright = yright + 2
            # ax.set_xlim(0.0, yright)
        else:
            yright = yright + 1
            # ax.set_xlim(0.0, yright)
        ax.set_xticks(np.arange(0.0, yright, 2.0))
        ax.set_xticklabels(np.arange(0.0, yright, 2.0, dtype=int))
    ## =====================================================================
    elif interval == 3:
        ## interval = 3
        if yright % 3 == 0:
            yright = yright + 1
            # ax.set_xlim(0.0, yright)
        else:
            yright = yright + 1
            while yright % 3 == 0:
                yright = yright + 1
                # ax.set_xlim(0.0, yright)
            # ax.set_xlim(0.0, yright)
        ax.set_xticks(np.arange(0.0, yright, 3.0))
        ax.set_xticklabels(np.arange(0.0, yright, 3.0, dtype=int))
    ## =========================================================================

    ax.set_xticks(np.arange(0.0, yright, 2.0))
    ax.set_xticklabels(np.arange(0.0, yright, 2.0, dtype=int))

    # loc = plticker.MultipleLocator(
    #     base=2.0)  # this locator puts ticks at regular intervals
    # ax.xaxis.set_major_locator(loc)

    ax.set_yscale('log')
    ax.grid(color='black', linestyle='--', linewidth=0.5)
    # ax.set_yticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1])
    ax.set_xlabel(r'Number of basis functions $\textit{N}$',
                  fontsize=font_size)
    ax.set_ylabel(r'Absolute error', fontsize=font_size)
    ax.tick_params(labelsize=tick_size)
    ax.legend(fontsize=font_size)
    # plt.show()
    ax.set_box_aspect(1)
    fig.savefig(fig_path, bbox_inches='tight', dpi=600)
    plt.close()


def plot_error_analysis_manual(
    current_directory,
    RB_folder,
    fig_path,
):
    plt.rcParams['text.usetex'] = True

    error_folder = os.path.join(current_directory, RB_folder,
                                'manual_error_analysis')
    error_data_list = os.listdir(error_folder)
    N_list = [i + 1 for i in range(len(error_data_list))]
    max_error_u_list, max_error_estimator_u_list = [], []
    gmean_error_u_list, gmean_error_estimator_u_list = [], []

    for plot_header in ["error", "error_estimator"]:
        for N in N_list:
            error_file_path = os.path.join(error_folder, '%s.csv' % N)
            error_data = pd.read_csv(error_file_path)
            mean_error = np.average(error_data[plot_header].values.tolist())
            max_error = np.amax(error_data[plot_header].values.tolist())
            if plot_header == "error":
                max_error_u_list.append(max_error)
                gmean_error_u_list.append(mean_error)
            elif plot_header == "error_estimator":
                max_error_estimator_u_list.append(max_error)
                gmean_error_estimator_u_list.append(mean_error)

    # original from here
    N = N_list

    fig, ax = plt.subplots(nrows=1,ncols=1)
    fig.set_figwidth(value_width)
    # 'max(error_estimator_u)' : max error bound
    ax.plot(N,
            max_error_estimator_u_list,
            ':bo',
            linewidth=3.0,
            markersize=6.5,
            label=r'Max $\Delta^{N}$')
    # 'max(error_u)' : max truth error
    ax.plot(N,
            max_error_u_list,
            '-bo',
            linewidth=1.0,
            markersize=8.5,
            label=r'Max $||e||_X$',
            mfc='none')
    # 'gmean(error_estimator_u)' : mean error bound
    ax.plot(N,
            gmean_error_estimator_u_list,
            ':r^',
            linewidth=3.0,
            markersize=6.5,
            label=r'Average $\Delta^{N}$')
    # 'gmean(error_u)' : mean truth error
    ax.plot(N,
            gmean_error_u_list,
            '-r^',
            linewidth=1.0,
            markersize=8.5,
            label=r'Average $||e||_X$',
            mfc='none')

    font_size = fontsize
    tick_size = fontsize
    ax.set_xticks(N)
    ax.set_yscale('log')
    ax.grid(color='black', linestyle='--', linewidth=0.5)
    # ax.set_yticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1])
    ax.set_xlabel(r'Number of basis functions $\textit{N}$',
                  fontsize=font_size)
    ax.set_ylabel(r'Absolute error', fontsize=font_size)
    ax.tick_params(labelsize=tick_size)
    xleft, yright = ax.get_xlim()
    interval = 3
    ## =====================================================================
    if interval == 2:
        ## interval = 2
        if yright % 2 == 0:
            yright = yright + 2
            # ax.set_xlim(0.0, yright)
        else:
            yright = yright + 1
            # ax.set_xlim(0.0, yright)
        ax.set_xticks(np.arange(0.0, yright, 2.0))
        ax.set_xticklabels(np.arange(0.0, yright, 2.0, dtype=int))
    ## =====================================================================
    elif interval == 3:
        ## interval = 3
        if yright % 3 == 0:
            yright = yright + 1
            # ax.set_xlim(0.0, yright)
        else:
            yright = yright + 1
            while yright % 3 == 0:
                yright = yright + 1
                # ax.set_xlim(0.0, yright)
            # ax.set_xlim(0.0, yright)
        ax.set_xticks(np.arange(0.0, yright, 3.0))
        ax.set_xticklabels(np.arange(0.0, yright, 3.0, dtype=int))
    ## =========================================================================


    ax.tick_params(labelsize=tick_size)
    # ax.legend(fontsize=16, framealpha=0.4, loc='upper right')
    ax.legend(ncol=2, bbox_to_anchor=(0.05, 1.2, 1., .102), fontsize=fontsize-3,
            fancybox=False, framealpha=0.8)
    ax.set_box_aspect(1)
    fig.savefig(fig_path, bbox_inches='tight', dpi=600)
    plt.close()


def plot_error_analysis_effectivity(RB_folder, solution_folder, fig_path):
    error_file_path = os.path.join(RB_folder, solution_folder,
                                   'error_analysis', 'error_analysis',
                                   'solution_u_error.csv')
    error_data = pd.read_csv(error_file_path, delimiter=';')
    N = error_data[['N']].values.tolist()
    max_error = error_data['max(error_u)'].values.tolist()
    max_bound = error_data['max(error_estimator_u)'].values.tolist()
    mean_error = error_data['gmean(error_u)'].values.tolist()
    mean_bound = error_data['gmean(error_estimator_u)'].values.tolist()
    eff_max = [max_bound[idx] / max_error[idx] for idx in range(len(N))]
    eff_avg = [mean_bound[idx] / mean_error[idx] for idx in range(len(N))]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    ax.plot(N,
            eff_max,
            '-bo',
            linewidth=1.0,
            markersize=8.5,
            label='Max $\eta^N(\mu)$',
            mfc='none')
    ax.plot(N,
            eff_avg,
            '-r^',
            linewidth=1.0,
            markersize=8.5,
            label='Average $\eta^N(\mu)$',
            mfc='none')
    font_size = fontsize
    tick_size = fontsize
    ax.set_xticks(N)
    # ax.set_yscale('log')
    ax.grid(color='black', linestyle='--', linewidth=0.5)
    # ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
    ax.set_xlabel(r'Number of basis functions $\textit{N}$',
                  fontsize=font_size)
    ax.set_ylabel(r'Effectivity factor $\eta^N(\mu)$', fontsize=font_size)
    ax.tick_params(labelsize=tick_size)
    ax.legend(fontsize=font_size, framealpha=0.4, loc='upper right')
    ax.set_box_aspect(1)
    fig.savefig(fig_path, bbox_inches='tight', dpi=600)
    plt.close()


# In[]

# def testing_SCM(problem, mu_range, folder_path, configuration):
#     print("Testing the SCM approach")
#     if os.path.isdir(folder_path):
#         pass
#     else:
#         os.makedirs(folder_path)

#     def N_generator(reduced_problem):
#         N = reduced_problem.N
#         if isinstance(N, dict):
#             N = min(N.values())
#         for n in range(1, N + 1):  # n = 1, ... N
#             yield n

#     saved_exact = []
#     saved_SCM = []
#     saved_E_f, saved_nu_f = [], []
#     saved_E_m, saved_nu_m = [], []

#     points = 5
#     E_f_list = np.linspace(min(mu_range[0]), max(mu_range[0]), points)
#     nu_f_list = np.linspace(min(mu_range[1]), max(mu_range[1]), points)
#     E_f_list = np.linspace(min(mu_range[2]), max(mu_range[2]), points)
#     nu_f_list = np.linspace(min(mu_range[3]), max(mu_range[3]), points)
#     i = 0
#     for nu_m in nu_f_list:
#         for E_m in E_f_list:
#             for nu_f in nu_f_list:
#                 for E_f in E_f_list:
#                     i += 1
#                     print(f"index: {i}")
#                     problem.set_mu((E_f, nu_f, E_m, nu_m))

#                     saved_SCM.append(
#                         problem.get_stability_factor_lower_bound())
#                     saved_exact.append(problem.evaluate_stability_factor())
#                     saved_E_f.append(E_f)
#                     saved_nu_f.append(nu_f)
#                     saved_E_m.append(E_m)
#                     saved_nu_m.append(nu_m)
#     data = {
#         "E_f": saved_E_f,
#         "nu_f": saved_nu_f,
#         "E_m": saved_E_m,
#         "nu_m": saved_nu_m,
#         "exact_alpha_LB": saved_exact,
#         "SCM_alpha_LB": saved_SCM
#     }
#     df = pd.DataFrame(data)
#     data_path = os.path.join(folder_path, f'SCM_{configuration}.csv')
#     df.to_csv(data_path, index=False)

# def manual_error_analysis(configuration, corrector, problem, reduction_method,
#                           reduced_problem, test_file_folder, folder_path):

#     def N_generator(reduced_problem):
#         N = reduced_problem.N
#         if isinstance(N, dict):
#             N = min(N.values())
#         for n in range(1, N + 1):  # n = 1, ... N
#             yield n

#     print(f"Manual Error Analysis: {configuration}-{corrector}")
#     if os.path.isdir(folder_path):
#         pass
#     else:
#         os.makedirs(folder_path)
#     reduction_method.testing_set.load(test_file_folder, "testing_set")
#     mu_set = reduction_method.testing_set
#     for n in N_generator(reduced_problem):
#         time_per_basis_start = time.time()
#         exact_list, bound_list, effectivity_list = [], [], []
#         exact_alpha_LB_list, SCM_alpha_LB_list = [], []
#         residual_norm_list, residual_norm_squared_list = [], []
#         E_f_list = [mu[0] for mu in mu_set]
#         nu_f_list = [mu[1] for mu in mu_set]
#         E_m_list = [mu[2] for mu in mu_set]
#         nu_m_list = [mu[3] for mu in mu_set]
#         time_online, time_alphalowerbound, time_exact, time_bound = [], [], [],[]
#         time_per_basis, time_per_mu = [], []

#         for mu in mu_set:
#             time_online_start = time.time()
#             reduced_problem.set_mu(mu)
#             reduced_problem.solve(n)
#             time_online.append(time.time() - time_online_start)

#             # Save the residual norm squared
#             residual_norm_squared = reduced_problem.get_residual_norm_squared()
#             residual_norm = sqrt(abs(residual_norm_squared))

#             residual_norm_squared_list.append(residual_norm_squared)
#             residual_norm_list.append(residual_norm)

#             # print(f'residual_norm_squared: {residual_norm_squared}')

#             # Save the alpha lower bound
#             time_alphalowerbound_start = time.time()
#             exact_alpha_LB = reduced_problem.evaluate_stability_factor()
#             # SCM_approximation.get_stability_factor_upper_bound()
#             SCM_alpha_LB = reduced_problem.get_stability_factor_lower_bound()
#             time_alphalowerbound.append(time.time() -
#                                         time_alphalowerbound_start)

#             # exact_alpha_LB = problem.evaluate_stability_factor(
#             # )  # exact alpha p.34 p44. 4.25
#             # SCM_alpha_LB = problem.get_stability_factor_lower_bound()
#             exact_alpha_LB_list.append(exact_alpha_LB)
#             SCM_alpha_LB_list.append(SCM_alpha_LB)

#             # Save the error bound and effectivity
#             time_exact_start = time.time()
#             exact = reduced_problem.compute_error()  # true (exact) error
#             time_exact.append(time.time() - time_exact_start)

#             time_bound_start = time.time()
#             bound = reduced_problem.estimate_error()  # error bound
#             time_bound.append(time.time() - time_bound_start)

#             effectivity = bound / exact
#             exact_list.append(exact)
#             bound_list.append(bound)
#             effectivity_list.append(effectivity)

#             time_per_mu.append(time.time() - time_online_start)

#         time_per_basis = (time.time() - time_per_basis_start)

#         data = {
#             "E_f_list": E_f_list,
#             "nu_f_list": nu_f_list,
#             "E_m_list": E_m_list,
#             "nu_m_list": nu_m_list,
#             "exact": exact_list,
#             "bound": bound_list,  # error_estimator
#             "effectivity": effectivity_list,
#             "exact_alpha_LB": exact_alpha_LB_list,
#             "SCM_alpha_LB": SCM_alpha_LB_list,
#             "residual_norm_squared": residual_norm_squared_list,
#             "residual_norm": residual_norm_list,
#             # "time_online": time_online,
#             # "time_alphalowerbound": time_alphalowerbound,
#             # "time_exact": time_exact,
#             # "time_error_bound": time_bound,
#             # "time_per_mu": time_per_mu,
#             # "time_per_basis": time_per_basis
#         }
#         df = pd.DataFrame(data)
#         data_path = os.path.join(folder_path, '%s.csv' % str(n))
#         df.to_csv(data_path, index=False)

# ==============================================================================
# location of the rb_problem source code without SCM
# /usr/local/lib/python3.8/dist-packages/RBniCS-0.1.dev1-py3.8.egg/rbnics/problems/elliptic/elliptic_rb_reduced_problem.py
# def manual_error_analysis(
#     corrector,
#     problem,
#     reduction_method,
#     reduced_problem,
#     test_file_folder,
#     folder_path,
#     sample_size_RB,
#     SCM=None,
# ):

#     def N_generator(reduced_problem):
#         N = reduced_problem.N
#         if isinstance(N, dict):
#             N = min(N.values())
#         for n in range(1, N + 1):  # n = 1, ... N
#             yield n

#     print(f"Manual Error Analysis: {corrector}")
#     if os.path.isdir(folder_path):
#         pass
#     else:
#         os.makedirs(folder_path)
#     reduction_method.testing_set.load(test_file_folder, "testing_set")
#     mu_set = reduction_method.testing_set

#     # print("mu_set:", mu_set)
#     dimension1 = np.shape(mu_set)[0]
#     r_list = [mu for mu in mu_set]
#     # print(f"r_list: {r_list}")
#     for n in N_generator(reduced_problem):
#         print(f"basis function {n}th")

#         time_per_basis_start = time.time()
#         # errorTruth_list, errorEstimator_list, effectivity_list = [], [], []
#         # errorTruth_alphaLB_list = []
#         # SCM_alphaLB_list = []
#         # residual_norm_squared_list = [], []

#         # errorTruth_list, errorEstimator_list, effectivity_list, errorTruth_alphaLB_list, \
#         #     SCM_alphaLB_list, residual_norm_squared_list = (
#         #     [] for i in range(6))
#         # time_online, time_alphaLB, time_errorTruth, \
#         #     time_bound, time_per_basis, time_per_mu = (
#         #     [] for i in range(6))
#         # if dimension1 == sample_size_RB:
#         #     r_list = [mu for mu in mu_set]
#         # else:
#         #     r_list = [mu[0] for mu in mu_set]
#         error_list, error_estimator_list, \
#             relative_error_list, relative_error_estimator_list, \
#             effectivity_list, residual_norm_squared_list, \
#             SCM_alphaLB_list, errorTruth_alphaLB_list = (
#             [] for i in range(8))

#         for mu in mu_set:
#             reduced_problem.set_mu(mu)
#             reduced_problem.solve(n)

#             error = reduced_problem.compute_error()
#             error_estimator = reduced_problem.estimate_error()
#             relative_error = reduced_problem.compute_relative_error()
#             relative_error_estimator = reduced_problem.estimate_relative_error(
#             )
#             effectivity = error_estimator / error
#             residual_norm_squared = reduced_problem.get_residual_norm_squared()

#             # Save truth error
#             error_list.append(error)
#             # Save error estimator
#             error_estimator_list.append(error_estimator)
#             # Save raltive error
#             relative_error_list.append(relative_error)
#             # Save relative error estimator
#             relative_error_estimator_list.append(relative_error_estimator)
#             # Save effectivity
#             effectivity_list.append(effectivity)
#             # Save residual norm squared
#             residual_norm_squared_list.append(residual_norm_squared)

#             # Save the alpha lower bound
#             if SCM != None:
#                 errorTruth_alphaLB = reduced_problem.evaluate_stability_factor(
#                 )
#                 # if SCM == None, the .get_stability_factor_lower_bound() will
#                 # simply return .evaluate_stability_factor()
#                 SCM_alphaLB = reduced_problem.get_stability_factor_lower_bound(
#                 )
#                 SCM_alphaLB_list.append(SCM_alphaLB)
#             else:
#                 errorTruth_alphaLB = reduced_problem.truth_problem.get_stability_factor_lower_bound(
#                 )
#             errorTruth_alphaLB_list.append(errorTruth_alphaLB)

#         data = {
#             "r_list": r_list,
#             "error": error_list,
#             "error_estimator": error_estimator_list,
#             "relative_error": relative_error,
#             "relative_error_estimator": relative_error_estimator,
#             "effectivity": effectivity_list,
#             "errorTruth_alphaLB": errorTruth_alphaLB_list,
#             # "SCM_alphaLB": SCM_alphaLB_list,
#             "residual_norm_squared": residual_norm_squared_list,
#         }
#         df = pd.DataFrame(data)
#         data_path = os.path.join(folder_path, '%s.csv' % str(n))
#         df.to_csv(data_path, index=False)

# 4params_wo_Radius, 4params_w_Radius, radiusOnly
    # ==============================================================================
    # ==============================================================================
def manual_error_analysis(
    configuration, corrector, 
    problem, reduction_method,
    reduced_problem, test_file_folder,folder_path,
    type_manual_error_analysis = "4params_w_Radius", 
    ):
    """
    # type_manual_error_analysis:=
        * 4params_wo_Radius 
        * 4params_w_Radius 
        * 5params_w_Radius 
        * radiusOnly
    """
    def N_generator(reduced_problem):
        N = reduced_problem.N
        if isinstance(N, dict):
            N = min(N.values())
        for n in range(1, N + 1):  # n = 1, ... N
            yield n
    
    if type_manual_error_analysis == "4params_wo_Radius": 

        print(f"Manual Error Analysis: {configuration}-{corrector}")
        if os.path.isdir(folder_path):
            pass
        else:
            os.makedirs(folder_path)
        reduction_method.testing_set.load(test_file_folder, "testing_set")
        mu_set = reduction_method.testing_set
        for n in N_generator(reduced_problem):
            time_per_basis_start = time.time()
            exact_list, bound_list, effectivity_list = [], [], []
            exact_alpha_LB_list, SCM_alpha_LB_list = [], []
            residual_norm_list, residual_norm_squared_list = [], []
            E_f_list = [mu[0] for mu in mu_set]
            nu_f_list = [mu[1] for mu in mu_set]
            E_m_list = [mu[2] for mu in mu_set]
            nu_m_list = [mu[3] for mu in mu_set]
            time_online, time_alphalowerbound, time_exact, time_bound = [], [], [],[]
            time_per_basis, time_per_mu = [], []

            for mu in mu_set:
                time_online_start = time.time()
                reduced_problem.set_mu(mu)
                reduced_problem.solve(n)
                time_online.append(time.time() - time_online_start)

                # Save the residual norm squared
                residual_norm_squared = reduced_problem.get_residual_norm_squared()
                residual_norm = sqrt(abs(residual_norm_squared))

                residual_norm_squared_list.append(residual_norm_squared)
                residual_norm_list.append(residual_norm)

                # print(f'residual_norm_squared: {residual_norm_squared}')

                # Save the alpha lower bound
                time_alphalowerbound_start = time.time()
                exact_alpha_LB = reduced_problem.evaluate_stability_factor()
                # SCM_approximation.get_stability_factor_upper_bound()
                SCM_alpha_LB = reduced_problem.get_stability_factor_lower_bound()
                time_alphalowerbound.append(time.time() -
                                            time_alphalowerbound_start)

                # exact_alpha_LB = problem.evaluate_stability_factor(
                # )  # exact alpha p.34 p44. 4.25
                # SCM_alpha_LB = problem.get_stability_factor_lower_bound()
                exact_alpha_LB_list.append(exact_alpha_LB)
                SCM_alpha_LB_list.append(SCM_alpha_LB)

                # Save the error bound and effectivity
                time_exact_start = time.time()
                exact = reduced_problem.compute_error()  # true (exact) error
                time_exact.append(time.time() - time_exact_start)

                time_bound_start = time.time()
                bound = reduced_problem.estimate_error()  # error bound
                time_bound.append(time.time() - time_bound_start)

                effectivity = bound / exact
                exact_list.append(exact)
                bound_list.append(bound)
                effectivity_list.append(effectivity)

                time_per_mu.append(time.time() - time_online_start)

            time_per_basis = (time.time() - time_per_basis_start)

            data = {
                "E_f_list": E_f_list,
                "nu_f_list": nu_f_list,
                "E_m_list": E_m_list,
                "nu_m_list": nu_m_list,
                "exact": exact_list,
                "bound": bound_list,  # error_estimator
                "effectivity": effectivity_list,
                "exact_alpha_LB": exact_alpha_LB_list,
                "SCM_alpha_LB": SCM_alpha_LB_list,
                "residual_norm_squared": residual_norm_squared_list,
                "residual_norm": residual_norm_list,
                # "time_online": time_online,
                # "time_alphalowerbound": time_alphalowerbound,
                # "time_exact": time_exact,
                # "time_error_bound": time_bound,
                # "time_per_mu": time_per_mu,
                # "time_per_basis": time_per_basis
            }
            df = pd.DataFrame(data)
            data_path = os.path.join(folder_path, '%s.csv' % str(n))
            df.to_csv(data_path, index=False)


    # ==============================================================================
    # 4 Parameters: ================================================================
    # (i) radius, (ii) ratio of Young modulus ======================================
    # (iii) fiber and (iv) matrix Poisson's ratio ==================================
    # ==============================================================================
    # def manual_error_analysis(
    #         configuration, corrector, 
    #         problem, reduction_method,
    #         reduced_problem, test_file_folder, 
    #         folder_path, list_range = 8):

    #     def N_generator(reduced_problem):
    #         N = reduced_problem.N
    #         if isinstance(N, dict):
    #             N = min(N.values())
    #         for n in range(1, N + 1):  # n = 1, ... N
    #             yield n

    elif type_manual_error_analysis ==  "4params_w_Radius":
        list_range = 8
        print(f"Manual Error Analysis: {configuration}-{corrector}")
        if os.path.isdir(folder_path):
            pass
        else:
            os.makedirs(folder_path)
        reduction_method.testing_set.load(test_file_folder, "testing_set")
        mu_set = reduction_method.testing_set
        for n in N_generator(reduced_problem):
            time_per_basis_start = time.time()
            exact_list, bound_list, effectivity_list = [], [], []
            exact_alpha_LB_list, SCM_alpha_LB_list = [], []
            residual_norm_list, residual_norm_squared_list = [], []
            # ======================================================================
            r_list = [mu[0] for mu in mu_set]
            ratio_E_list = [mu[1] for mu in mu_set]
            nu_f_list = [mu[2] for mu in mu_set]
            nu_m_list = [mu[3] for mu in mu_set]
            # ======================================================================
            time_online, time_alphalowerbound, time_exact, time_bound = [], [], [],[]
            time_per_basis, time_per_mu = [], []

            error_list, error_estimator_list, \
                relative_error_list, relative_error_estimator_list, \
                effectivity_list, residual_norm_squared_list, \
                SCM_alphaLB_list, errorTruth_alphaLB_list = (
                [] for i in range(list_range))

            for mu in mu_set:
                time_online_start = time.time()
                reduced_problem.set_mu(mu)
                reduced_problem.solve(n)
                time_online.append(time.time() - time_online_start)

                # Save the residual norm squared
                residual_norm_squared = reduced_problem.get_residual_norm_squared()
                residual_norm = sqrt(abs(residual_norm_squared))

                residual_norm_squared_list.append(residual_norm_squared)
                residual_norm_list.append(residual_norm)

                # print(f'residual_norm_squared: {residual_norm_squared}')

                # Save the alpha lower bound
                time_alphalowerbound_start = time.time()
                exact_alpha_LB = reduced_problem.evaluate_stability_factor()
                # SCM_approximation.get_stability_factor_upper_bound()
                SCM_alpha_LB = reduced_problem.get_stability_factor_lower_bound()
                time_alphalowerbound.append(time.time() -
                                            time_alphalowerbound_start)

                # exact_alpha_LB = problem.evaluate_stability_factor(
                # )  # exact alpha p.34 p44. 4.25
                # SCM_alpha_LB = problem.get_stability_factor_lower_bound()
                exact_alpha_LB_list.append(exact_alpha_LB)
                SCM_alpha_LB_list.append(SCM_alpha_LB)

                # Save the error bound and effectivity
                time_exact_start = time.time()
                exact = reduced_problem.compute_error()  # true (exact) error
                time_exact.append(time.time() - time_exact_start)

                time_bound_start = time.time()
                bound = reduced_problem.estimate_error()  # error bound
                time_bound.append(time.time() - time_bound_start)

                effectivity = bound / exact
                exact_list.append(exact)
                bound_list.append(bound)
                effectivity_list.append(effectivity)

                time_per_mu.append(time.time() - time_online_start)

            time_per_basis = (time.time() - time_per_basis_start)

            data = {
                "r_list": r_list,
                "ratio_E_list": ratio_E_list,
                "nu_f_list": nu_f_list,
                "nu_m_list": nu_m_list,
                "exact": exact_list,
                "bound": bound_list,  # error_estimator
                "effectivity": effectivity_list,
                "exact_alpha_LB": exact_alpha_LB_list,
                "SCM_alpha_LB": SCM_alpha_LB_list,
                "residual_norm_squared": residual_norm_squared_list,
                "residual_norm": residual_norm_list,
                "time_online": time_online,
                # "time_alphalowerbound": time_alphalowerbound,
                "time_exact": time_exact,
                "time_error_bound": time_bound,
                # "time_per_mu": time_per_mu,
                # "time_per_basis": time_per_basis
            }
            # data[time_online] = time_online
            df = pd.DataFrame(data)
            data_path = os.path.join(folder_path, '%s.csv' % str(n))
            df.to_csv(data_path, index=False)

    # ==============================================================================
    # ==============================================================================
    # ==============================================================================
    elif type_manual_error_analysis == "radiusOnly":
        SCM=None

        print(f"Manual Error Analysis: {corrector}")
        if os.path.isdir(folder_path):
            pass
        else:
            os.makedirs(folder_path)
        reduction_method.testing_set.load(test_file_folder, "testing_set")
        mu_set = reduction_method.testing_set

        # print("mu_set:", mu_set)
        dimension1 = np.shape(mu_set)[0]
        r_list = [mu for mu in mu_set]
        # print(f"r_list: {r_list}")
        for n in N_generator(reduced_problem):
            print(f"basis function {n}th")

            time_per_basis_start = time.time()
            # errorTruth_list, errorEstimator_list, effectivity_list = [], [], []
            # errorTruth_alphaLB_list = []
            # SCM_alphaLB_list = []
            # residual_norm_squared_list = [], []

            # errorTruth_list, errorEstimator_list, effectivity_list, errorTruth_alphaLB_list, \
            #     SCM_alphaLB_list, residual_norm_squared_list = (
            #     [] for i in range(6))
            # time_online, time_alphaLB, time_errorTruth, \
            #     time_bound, time_per_basis, time_per_mu = (
            #     [] for i in range(6))
            # if dimension1 == sample_size_RB:
            #     r_list = [mu for mu in mu_set]
            # else:
            #     r_list = [mu[0] for mu in mu_set]
            error_list, error_estimator_list, \
                relative_error_list, relative_error_estimator_list, \
                effectivity_list, residual_norm_squared_list, \
                SCM_alphaLB_list, errorTruth_alphaLB_list = (
                [] for i in range(8))

            data = {}

            for mu in mu_set:
                # print(f"mu: {mu}")
                reduced_problem.set_mu(mu)
                reduced_problem.solve(n)

                error = reduced_problem.compute_error()
                error_estimator = reduced_problem.estimate_error()
                relative_error = reduced_problem.compute_relative_error()
                relative_error_estimator = reduced_problem.estimate_relative_error()
                effectivity = error_estimator / error
                residual_norm_squared = reduced_problem.get_residual_norm_squared()

                # Save truth error
                error_list.append(error)
                # Save error estimator
                error_estimator_list.append(error_estimator)
                # Save raltive error
                relative_error_list.append(relative_error)
                # Save relative error estimator
                relative_error_estimator_list.append(relative_error_estimator)
                # Save effectivity
                effectivity_list.append(effectivity)
                # Save residual norm squared
                residual_norm_squared_list.append(residual_norm_squared)

                # Save the alpha lower bound
                if SCM != None:
                    errorTruth_alphaLB = reduced_problem.evaluate_stability_factor()
                    # if SCM == None, the .get_stability_factor_lower_bound() will
                    # simply return .evaluate_stability_factor()
                    SCM_alphaLB = reduced_problem.get_stability_factor_lower_bound()
                    SCM_alphaLB_list.append(SCM_alphaLB)
                    data["SCM_alphaLB"]= SCM_alphaLB_list
                else:
                    errorTruth_alphaLB = reduced_problem.truth_problem.get_stability_factor_lower_bound()
                errorTruth_alphaLB_list.append(errorTruth_alphaLB)

            data["r_list"] =  r_list
            data["error"] = error_list
            data["error_estimator"] = error_estimator_list
            data["relative_error"] = relative_error
            data["relative_error_estimator"] = relative_error_estimator
            data["effectivity"] = effectivity_list
            data["errorTruth_alphaLB"] = errorTruth_alphaLB_list
            data["residual_norm_squared"] = residual_norm_squared_list
            # data = {
            #     "r_list": r_list,
            #     "error": error_list,
            #     "error_estimator": error_estimator_list,
            #     "relative_error": relative_error,
            #     "relative_error_estimator": relative_error_estimator,
            #     "effectivity": effectivity_list,
            #     "errorTruth_alphaLB": errorTruth_alphaLB_list,
            #     # "SCM_alphaLB": SCM_alphaLB_list,
            #     "residual_norm_squared": residual_norm_squared_list,
            # }
            df = pd.DataFrame(data)
            data_path = os.path.join(folder_path, '%s.csv' % str(n))
            df.to_csv(data_path, index=False)
# ==============================================================================
# ==============================================================================


class error_estimation():
    # """# 5. posteriori error estimation """

    def __init__(
        self,
        corrector,
        modelName,
        mu_range,
        problem,
        sample_size_RB,
        reduction_method,
        reduced_problem,
        EIM_testing_set=None,
        **kwargs,
    ):
        self.corrector = corrector
        self.modelName = modelName
        self.sample_size_RB = sample_size_RB
        self.EIM_testing_set = EIM_testing_set

        self.problem = problem
        self.reduction_method = reduction_method
        self.reduced_problem = reduced_problem
        self.mu_range = mu_range
        # posteriori_error_start = time.time()

        self.test_folder_path = os.path.join(os.getcwd(), f"{modelName}")
        self.test_file_folder = os.path.join(self.test_folder_path,
                                             "testing_set")
        self.folder_path = os.path.join(self.test_folder_path,
                                        "manual_error_analysis")

        # 5.1 Perform a manual error analysis
        print(f"5.1 Perform a manual error analysis: {modelName}")

        # self.reduction_method.initialize_testing_set(
        #     self.sample_size_RB,
        #     DEIM=self.EIM_testing_set,
        #     enable_import=True,
        # )
        # print("Latin Hypercube sampling for testing_set")
        # self.sample_size_RB = Latin_Hypercube_sampling_func(
        #     sample_size_RB=self.sample_size_RB,
        #     modelName=modelName,
        #     mu_range=self.mu_range,
        #     type_of_sample_set="testing_set",
        # )

        # exec(open("./LatinHyperCube_sampling.py").read())

        # ----------------------------------------------------------------------
        # 5.1 Perform a manual error analysis
        self._manualErrorAnalysis()
        # ----------------------------------------------------------------------
        # Perform an built-in error analysis to estimate
        # self.callErrorAnalysis()

    def _manualErrorAnalysis(self, ):
        # time_manual_error_analysis_start = time.time()
        manual_error_analysis(
            corrector=self.corrector,
            problem=self.problem,
            reduced_problem=self.reduced_problem,
            reduction_method=self.reduction_method,
            test_file_folder=self.test_folder_path,
            folder_path=self.folder_path,
            sample_size_RB=self.sample_size_RB,
        )
        # time_series["manual_error_analysis_{corrector}"] = (
        #     time.time() - time_manual_error_analysis_start)

        # # posteriori_error_end = time.time()
        # print(
        #     f"- error analysis - {configuration} - {corrector}: {time.time() - posteriori_error_start} s-"
        # )
        # time_series["total error analysis {corrector}"] = (
        #     time.time() - time_manual_error_analysis_start)

    def callErrorAnalysis(self, ):
        # 5.2 Perform an built-in error analysis to estimate
        # print("5.2 Perform an built-in error analysis to estimate")
        self.reduction_method.error_analysis(filename=f"error_analysis")

# ==============================================================================
# ==============================================================================




# In[]

# ==============================================================================
# ==============================================================================
def basis_function_grad(V, reduced_problem, grad, modelName):
    _basis_function_grad = copy.copy(reduced_problem.basis_functions)
    print(reduced_problem.N)
    for i in range(reduced_problem.N):
        """
        input V are vector-space
        _basis_function_grad are
          * gradient of basis_i over dx: 
          where 
            * grad = 0 => dx1 
            * grad = 1 => dx2
          * a vector-valued
        """
        basis_i = reduced_problem.basis_functions[i]
        # _basis_i_grad_dx1 = project(basis_0.dx(0), V)

        if grad == 0:  #  gradient of basis_i over dx1
            _basis_function_grad[i] = project(basis_i.dx(0), V)
            # basis_function_grad_0[i] = project(grad(basis_i)[0], V)
        if grad == 1:  #  gradient of basis_i over dx1
            _basis_function_grad[i] = project(basis_i.dx(1), V)
        # ======================================================================

    # if False:  # this is testing are - Test work on 17:38, Apr 9
    #     # ==========================================================================
    #     # calculate the vector-valued gradient over dx 1 or 2 (based on input grad)
    #     grad_u_dx = _basis_function_grad * reduced_problem._solution

    #     # ==========================================================================
    #     print("plot grad_u_dx")
    #     plt.figure()
    #     plt.jet()
    #     plt.colorbar(
    #         plot(grad_u_dx[0],
    #              reduced_problem=reduced_problem,
    #              title=f"grad_u_1_dx{grad+1}"))
    #     plt.show()
    #     plt.close()

    #     # ==========================================================================
    #     print("plot grad_u_dx")
    #     plt.figure()
    #     plt.jet()
    #     plt.colorbar(
    #         plot(grad_u_dx[1],
    #              reduced_problem=reduced_problem,
    #              title=f"grad_u_1_dx{grad+1}"))
    #     plt.show()
    #     plt.close()
    #     # u_sol = reduced_problem_11._solution
    #     # solution_RB = reduced_problem_11.basis_functions[:reduced_problem_11._solution.N] * reduced_problem_11._solution
    #     # ==========================================================================

    #     # ==========================================================================
    """ save - load reduced basis functions"""
    directory = os.path.join(os.getcwd(), modelName, f"basis_grad_dx{grad+1}")
    check_create_dir(directory)
    ## this worked
    reduced_problem.basis_functions.save(directory=directory, filename="basis")

    return _basis_function_grad