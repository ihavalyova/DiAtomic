import os
import random
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline
# from collections import defaultdict
from utils import Utils
from constants import Const
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

import matplotlib
matplotlib.use('TkAgg')


class Plotting:

    def __init__(self):
        pass

    @classmethod
    def _define_main_keywords(cls):

        return {
            0: 'what',
            1: 'mlevels',
            2: 'fformat',
            3: 'show',
            4: 'path',
            5: 'props',
        }

    @classmethod
    def _define_what_keywords(cls):

        return {
            'residuals_level': cls.plot_residuals_level,
            'residuals': cls.full_residuals_plot,
            'potentials_on_grid': cls.plot_potentials_on_grid,
            'couplings_on_grid': cls.plot_couplings_on_grid,
        }

    @staticmethod
    def plot(cls, **kwargs):

        keys = cls._define_main_keywords()
        what_keys = cls._define_what_keywords()

        what = kwargs[keys[0]].lower()
        mlevels = kwargs.get(keys[1])
        frm = kwargs.get(keys[2]) or 'png'
        show = kwargs.get(keys[3]) or False
        path = kwargs.get(keys[4])
        props = kwargs.get(keys[5])

        if what not in what_keys.keys():
            raise SystemExit(f'Error: Invalid plotting parameter {what}.')

        what_keys[what](
            mlevels=mlevels, path=path, fformat=frm, show=show, props=props
        )

    @classmethod
    def full_residuals_plot(cls, mlevels, path=None, fformat='png',
                            show=False, props=None):

        path = path or Utils.get_plot_dir('resid_path')

        data = mlevels.out_data

        gen_color = cls._get_color()
        gen_marker = cls._get_marker(limit=True)

        fillstyle = ['none', 'none', 'full', 'full']
        random.shuffle(fillstyle)

        props_def = {
            'xlabel': 'Energy',
            'ylabel': 'E_obs - E_calc',
            'title': 'Residuals',
        }

        # will change props values
        props = cls._combine_props(props, props_def)

        nisotopes = len(mlevels.masses)

        for _ in range(1, nisotopes+1):
            mrng = (9*(nisotopes-1)+1, (nisotopes+9)-1)

            # split e- and f-levels
            if mrng is None:
                edata = data[data[:, 6] == 1]
                fdata = data[data[:, 6] == 0]
            else:
                edata = data[
                    data[:, 6] == 1 & (mrng[0] <= data[:, 7].all() <= mrng[1])
                ]
                fdata = data[
                    data[:, 6] == 0 & (mrng[0] <= data[:, 7].all() <= mrng[1])
                ]

            line_styles = {
                'color': (next(gen_color), next(gen_color)),
                'marker': (next(gen_marker), next(gen_marker)),
                'markersize': (6.5, 6.5),
                'markeredgewidth': (0.8, 0.8),
                'linewidth': (0.0, 0.0),
                'fillstyle': ('none', 'none')
            }

            fig, ax = plt.subplots()

            ax.plot(
                edata[:, 9],
                edata[:, 10],
                color=line_styles['color'][0],
                marker=line_styles['marker'][0],
                markersize=line_styles['markersize'][0],
                markeredgewidth=line_styles['markeredgewidth'][0],
                linewidth=line_styles['linewidth'][0],
                fillstyle=line_styles['fillstyle'][0]
            )

            ax.plot(
                fdata[:, 9],
                fdata[:, 10],
                color=line_styles['color'][1],
                marker=line_styles['marker'][1],
                markersize=line_styles['markersize'][1],
                markeredgewidth=line_styles['markeredgewidth'][1],
                linewidth=line_styles['linewidth'][1],
                fillstyle=line_styles['fillstyle'][1]
            )

            # annotate average uncertanty
            average_uncert = np.sum(data[:, 11]) / data[:, 11].shape

            plt.axhline(
                y=average_uncert, color='r', linestyle='-', linewidth=0.5
            )
            plt.axhline(
                y=-1.0*average_uncert, color='r', linestyle='-', linewidth=0.5
            )

            os.makedirs(path, exist_ok=True)

            ax.set_title(props['title'])
            ax.set_xlabel(props['xlabel'])
            ax.set_ylabel(props['ylabel'])
            ax.legend(['e-parity levels', 'f-parity levels'], loc=0)
            # fig.tight_layout()

            fig.set_size_inches(13, 7)
            fig_path = os.path.join(path, 'residuals' + f'.{fformat}')
            fig.savefig(fig_path, format=fformat, dpi=100)  # transparent=True
            print(f'Created figure: {fig_path}', sep='')
            if show:
                plt.show()

    @classmethod
    def plot_couplings_on_grid(cls, mlevels, path=None, fformat='png',
                               show=False, props=None):

        path = path or Utils.get_plot_dir('func_path')

        props_def = {
            'xlabel': 'Internuclear Distance',
            'ylabel': 'Interaction',
            'title': 'Couplings',
        }

        # will change props values
        props = cls._combine_props(props, props_def)

        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        fig4, ax4 = plt.subplots()

        gen_color = cls._get_color()
        gen_marker = cls._get_marker()

        fgrid_cols = np.hstack((
            mlevels.rgrid[:, np.newaxis],
            mlevels.fgrid.reshape(mlevels.ncp, mlevels.ngrid).T
        ))
        rgrid_col = fgrid_cols[:, 0] * Const.bohr
        n = rgrid_col.shape[0]

        leg1, leg2, leg3, leg4 = [], [], [], []

        for col in range(1, fgrid_cols.shape[1]):
            markers_on = [0, int(n/4), int(n/2), int(n/1.6), int(n/1.2), -1]

            if 'spin-orbit' in mlevels.couplings[col-1].coupling:
                ax1.plot(
                    rgrid_col,
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )
                ax1.set_title('Spin-Orbit interactions')
                ax1.set_xlabel(props['xlabel'])
                ax1.set_ylabel(props['ylabel'])
                ax1.set_ylim(-800, -200)  # only for nih - to remove
                leg1.append(mlevels.couplings[col-1].interact)

            if 'LJ' in mlevels.couplings[col-1].coupling:
                ax2.plot(
                    rgrid_col,
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )

                ax2.set_title('LJ interactions')
                ax2.set_xlabel(props['xlabel'])
                ax2.set_ylabel(props['ylabel'])
                ax2.set_ylim(0, 2)
                leg2.append(mlevels.couplings[col-1].interact)

            if 'SJ' in mlevels.couplings[col-1].coupling:
                ax3.plot(
                    rgrid_col,
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )

                ax3.set_title('SJ interactions')
                ax3.set_xlabel(props['xlabel'])
                ax3.set_ylabel(props['ylabel'])
                ax3.set_ylim(0, 3.5)
                leg3.append(mlevels.couplings[col-1].interact)

            if 'LambdaD' in mlevels.couplings[col-1].coupling:
                ax4.plot(
                    rgrid_col,
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )

                ax4.set_title('Lambda doubling interactions')
                ax4.set_xlabel(props['xlabel'])
                ax4.set_ylabel(props['ylabel'])
                ax4.set_ylim(-0.1, 0.1)
                leg4.append(mlevels.couplings[col-1].interact)

        ax1.legend(leg1, loc=0)
        ax2.legend(leg2, loc=0)
        ax3.legend(leg3, loc=0)
        ax4.legend(leg4, loc=0)

        os.makedirs(path, exist_ok=True)

        fig1_path = os.path.join(path, f'SO.{fformat}')
        fig2_path = os.path.join(path, f'LJ.{fformat}')
        fig3_path = os.path.join(path, f'SJ.{fformat}')
        fig4_path = os.path.join(path, f'LD.{fformat}')

        cls.save_figure(fig1, fig1_path, fformat, "tight")
        cls.save_figure(fig2, fig2_path, fformat, "tight")
        cls.save_figure(fig3, fig3_path, fformat, "tight")
        cls.save_figure(fig4, fig4_path, fformat, "tight")

        if show:
            plt.show()

    @classmethod
    def save_figure(self, fig, path, fformat, bbox):

        fig.savefig(path, format=fformat, bbox_inches=bbox)

    @classmethod
    def plot_residuals_level(cls, mlevels, path=None, fformat='png',
                             show=False, props=None):

        plt.rcParams.update({'font.size': 8})
        res_path = Utils.get_plot_dir('res_path')
        os.makedirs(res_path, exist_ok=True)
        path = path or res_path

        props_def = {
            'xlabel': 'J',
            'ylabel': 'E_obs - E_calc, cm-1',
            'title': 'v={}, state={}, isotop={}',
        }

        # will change props values
        props = cls._combine_props(props, props_def)
        data = mlevels.out_data

        # split data by v, then split by state, then extract
        # marker ranges and get the data for each isotope
        vsplit = np.split(data, np.where(np.diff(data[:, 1]))[0]+1)

        # TODO: find more clever way
        plt.rcParams.update({'figure.max_open_warning': 0})

        for v, _ in enumerate(vsplit):

            vn = int(vsplit[v][:, 1][0])

            vsplit[v] = vsplit[v][vsplit[v][:, -1].argsort()]

            ssplit = np.split(
                vsplit[v], np.where(np.diff(vsplit[v][:, -1]))[0]+1
            )

            for sti, _ in enumerate(ssplit):

                stn = int(ssplit[sti][:, -1][0])

                mrange = cls._map_marker(ssplit[sti][:, 7])

                for i, m in enumerate(mrange):
                    splitted_data = ssplit[sti][m, :]

                    fdata = splitted_data[splitted_data[:, 6] == 0]
                    edata = splitted_data[splitted_data[:, 6] == 1]

                    fig, ax_data = plt.subplots()

                    ax_data.plot(
                        edata[:, 2],
                        -edata[:, 10],
                        color='darkblue',
                        marker='o',
                        markersize=6,
                        markeredgewidth=0.7,
                        fillstyle='none',
                        linewidth=0,
                    )

                    ax_data.plot(
                        fdata[:, 2],
                        -fdata[:, 10],
                        color='darkgreen',
                        marker='X',
                        markersize=6,
                        markeredgewidth=0.7,
                        fillstyle='none',
                        linewidth=0,
                    )

                    ax_data.set_xlabel(props['xlabel'])
                    ax_data.set_ylabel(props['ylabel'])
                    ax_data.legend(['e', 'f'], loc=0, prop={'size': 6})
                    ax_data.set_title(props['title'].format(vn, stn, i+1))

                    ax_data.xaxis.set_major_locator(MultipleLocator(20))
                    ax_data.xaxis.set_major_formatter(FormatStrFormatter('%d'))
                    ax_data.xaxis.set_minor_locator(MultipleLocator(5))
                    ax_data.grid(which='both')
                    ax_data.grid(which='major', alpha=0.5)
                    ax_data.grid(which='minor', alpha=0.2)

                    fig.set_size_inches(7, 5)

                    figname = f'residual_{vn}_{stn}_{i}.{fformat}'
                    figpath = os.path.join(path, figname)
                    fig.savefig(
                        figpath, format=fformat, dpi=100, transparent=False
                    )
                    fig.tight_layout()

                    print(f'Created figure: {figpath}', sep='')

                    if show:
                        plt.show()

                    # plt.cla()
                    # plt.gcf().clear()
        # plt.close(fig)

    @classmethod
    def plot_potentials_on_grid(cls, mlevels, path=None, fformat='png',
                                show=False, props=None):

        path = path or Utils.get_plot_dir('plot_path')

        # to add x lim and y lim
        props_def = {
            'xlabel': 'Internuclear Distance',
            'ylabel': 'Potential energy',
            'title': 'Potentials',
        }

        # will change props values
        props = cls._combine_props(props, props_def)

        _, ax = plt.subplots()

        unique_pfiles = set()
        for channel in mlevels.channels:
            unique_pfiles.add(channel.filep)

        # need to keep the unique file names in ordered data structue
        unique_pfiles = list(unique_pfiles)

        gen_color = cls._get_color()

        for i, _ in enumerate(unique_pfiles):
            gridr = mlevels.rgrid * Const.bohr
            gridu_range = mlevels.ugrid[i*mlevels.ngrid:(i+1)*mlevels.ngrid]
            gridu = gridu_range * Const.hartree

            ax.plot(
                gridr,
                gridu,
                color=next(gen_color),
                markersize=0,
                markeredgewidth=0.7,
                linewidth=1.5,
            )

        ax.set_title(props['title'])
        ax.set_xlabel(props['xlabel'])
        ax.set_ylabel(props['ylabel'])
        ax.legend(unique_pfiles, loc=0)

        figpath = os.path.join(path, f'potentials.{fformat}')
        plt.savefig(figpath, format=fformat)

        print(f'Created figure: {figpath}', sep='')

        if show:
            plt.show()

    @classmethod
    def plot_potentials_points(cls, files, path=None, fformat='png',
                               show=False, ipoints=None, xlim=None, ylim=None):
        # cannot not be called as option from plot()
        # only for plotting pointwise potentials

        _, ax = plt.subplots()

        path = path or Utils.get_plot_dir('plot_path')

        gen_color = cls._get_color()
        gen_marker = cls._get_marker()
        xlim = xlim or (None, None)
        ylim = ylim or (None, None)

        is_failed = True

        for pfile in files:
            try:
                pdata = np.loadtxt(pfile, skiprows=1)
                x, y = pdata[:, 0], pdata[:, 1]

                ipoints = ipoints or x.shape[0]
                x_interp = np.linspace(x[0], x[-1], ipoints, endpoint=True)

                cs = CubicSpline(x, y, bc_type='natural')
                y_interp = cs(x_interp)

                ax.plot(
                    x_interp,
                    y_interp,
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=5.5,
                    markeredgewidth=0.4,
                    linewidth=1.5,
                    fillstyle='none'
                )
                # should pass at least once
                is_failed = False
            except Exception as err:
                print(
                    f'Error: The file {pfile} does not exist or is not in '
                    f'the correct format and cannot be poltted!\n{err}'
                )

        if not is_failed:

            ax.set_title('Potentials')
            ax.set_xlabel('Internuclear distance')
            ax.set_ylabel('Potential energy')
            ax.set_xlim(xlim[0], xlim[1])
            ax.set_ylim(ylim[0], ylim[1])
            ax.legend(files, loc=0)

            figpath = os.path.join(path, f'potential_points.{fformat}')
            plt.savefig(figpath, format=fformat)

            print(f'Created figure: {figpath}', sep='')

            if show:
                plt.show()

    @staticmethod
    def hcolormesh(mlevels, rows=None, cols=None, path=None,
                   fformat='png', show=False):

        """create a colormesh of the Hamiltonian matrix

        Args:
            mlevels (object): MoleculeLevels object
            rows (tuple, optional): the row indecies. Defaults to None.
            cols (tuple, optional): the col indices. Defaults to None.
            path (str, optional): path where to save the created figure.
            Defaults to None.
            fformat (str, optional): File extension. Defaults to 'png'.
            show (bool, optional): whether to show the image.
            Defaults to False.

        Remarks:
            1. rows and cols should be tuples with two integer numbers -
            the first and the last row and the first and the last column
            of the submatrix whcih to be plotted
            2. cannot be called as option from plot()
        """

        # TODO: The Hamiltonian is upper diagonal matrix - convert it to full
        hmat = np.copy(mlevels.hmatrix)
        hmat[hmat == 0.0] = np.nan

        rows = rows or (0, hmat.shape[0])
        cols = cols or (0, hmat.shape[1])

        fig, ax = plt.subplots()

        # im = ax.pcolormesh(hmat[rows[0]:rows[1], cols[0]:cols[1]],
        #   vmin=1.0e-8, vmax=0.01, cmap='RdBu_r')

        '''
        vmin and vmax set the normalization range.
        By default scale scalar data to the [0, 1] range
        '''

        im = ax.matshow(
            hmat[rows[0]:rows[1], cols[0]:cols[1]], vmin=1.0e-8, vmax=0.01,
            cmap='RdBu_r', aspect='auto', interpolation=None
        )

        fig.colorbar(im)
        ax.set_title(
            f'Hamiltonian matrix colormesh: '
            f'rows={rows[0]}:{rows[1]}, cols={cols[0]}:{cols[1]}'
        )

        # fig.tight_layout()
        path = path or Utils.get_plot_dir('plot_path')

        figpath = os.path.join(path, f'hcolormesh.{fformat}')
        plt.savefig(figpath, format=fformat)

        print(f'Created figure: {figpath}', sep='')

        if show:
            plt.show()

    @staticmethod
    def plot_wavefunctions(wffile, props=None, path=None,
                           subplots=True, fformat='png', title='',
                           xlabel='R', ylabel='wavefunction', leg=''):
        """
        Plot the interpolated wavefunction over the interpolation grid

        Args:
            wffile: file where the interpolated wavefunction is stored
            props: various properties concerning the style of the figure
            path: the directory path where the figure will be saved
            fformat: the file format of the figure (e.g. png, eps,...)
            title: the title of the figure
            subplots: whether the wavefunctions to be plotted as subplots
        Returns: None
        """

        wfdata = np.loadtxt(wffile)
        igrid = wfdata[:, 0]
        wavefunc = wfdata[:, 0:]

        path = path or Utils.get_plot_dir('wavefunc_path')

        os.makedirs(path, exist_ok=True)

        if not subplots:

            fig_name = os.path.join(path, 'wavefunctions' + f'.{fformat}')

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend(leg)
            plt.plot(igrid, wavefunc)
            plt.savefig(fig_name, format=fformat, dpi=250)

            print(f'Created figure: {fig_name}', sep='')
        else:
            # fig, ax = plt.subplots()
            # fig = plt.figure()
            # colors = self._get_color()

            for i in range(1, wavefunc.shape[1]+1):
                ax = plt.subplot(2, 3, i)
                ax.set_title(f'v={i-1}')
                ax.set_ylabel(ylabel)
                ax.set_xlabel(xlabel)
                plt.axis('on')

                plt.plot(igrid, wavefunc[:, i-1])

            plt.subplots_adjust(wspace=0.4, hspace=0.4)

            plt.show()

    @classmethod
    def _combine_props(cls, props, props_def):

        if props is None:
            return props_def
        else:
            for key, value in props_def.items():
                try:
                    props[key]
                except KeyError:
                    props[key] = value
            return props

    @classmethod
    def _map_marker(cls, markers):

        ranges = []
        for i in range(0, 60, 10):  # up to 6 isotopes
            res = np.where(np.logical_and(markers >= i, markers <= i+10))[0]
            if res.shape[0] > 0:
                ranges.append(res)

        return ranges

    @classmethod
    def _get_color(cls):

        bcolors = [
            '#1C03FF', '#4132C7', '#130B5B', '#10286B',
            '#355098', '#380AA2', '#3061B5', '#1C8FC9',
            '#00ABFF', '#10C7EC', '#10ECEC', '#BA5FFF',
            '#800080', '#1E90FF', '#4682B4', '#4B0082',
        ]
        random.shuffle(bcolors)

        rcolors = [
            '#FA0A0A', '#D44141', '#B53232', '#EF262D',
            '#D50D50', '#FC697C', '#FC69E3', '#A11589',
            '#A8639C', '#DF211A', '#F3592A', '#FF8F44',
            '#C71585', '#D2691E', '#DC143C', '#800000'
        ]
        random.shuffle(rcolors)

        gcolors = [
            '#00FF2B', '#3AA44C', '#165F22', '#54FC70',
            '#3EE22C', '#15730B', '#4EA644', '#65EE16',
            '#13987D', '#32B30B', '#77D11D', '#D1BE15',
            '#556B2F', '#2E8B57', '#008B8B', '#2E8B57'
        ]
        random.shuffle(gcolors)

        colors = []
        for color in zip(gcolors, bcolors, rcolors):
            g, b, r = color
            colors.append(g)
            colors.append(b)
            colors.append(r)

        for color in colors:
            yield color

    @classmethod
    def _get_marker(cls, limit=False):

        markers = [
            'o', 'v', '<', '>', 's', 'p', '*',
            '+', 'x', 'D', 'd', 'X', '^', 'h',
            '1', '2', '3', '4', 'o', 'v', 'X',
            'x', 's', '<', '>', 'D', 'd', '*',
            'o', 'v', 'X', 'x'
        ]

        if limit:
            markers = [
                'o', 'v', '<', '>', 's',
                '*', 'x', 'D', 'd', 'X'
            ]

        random.shuffle(markers)

        for marker in markers:
            yield marker
