import os
import random
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline
# from collections import defaultdict
from utils import Utils
from constants import Const

import matplotlib
matplotlib.use('TkAgg')


class Plotting:

    def __init__(self):
        self.plot_dir = 'plotting'
        self.func_dir = 'interactions'
        self.resid_dir = 'full_residuals'
        self.wavefunc_dir = 'wavefunctions'
        self.res_dir = 'residuals'

        self.plot_path = os.path.join(Utils.get_current_dir(), self.plot_dir)
        self.func_path = os.path.join(self.plot_path, self.func_dir)
        self.resid_path = os.path.join(self.plot_path, self.resid_dir)
        self.wavefunc_path = os.path.join(self.plot_path, self.wavefunc_dir)
        self.res_path = os.path.join(self.plot_path, self.res_dir)

        os.makedirs(self.plot_path, exist_ok=True)

    def _define_main_keywords(self):

        return {
            0: 'what',
            1: 'mlevels',
            2: 'fformat',
            3: 'show',
            4: 'path',
            5: 'props',
        }

    def _define_what_keywords(self):

        return {
            'residuals_level': self.plot_residuals_level,
            'residuals': self.full_residuals_plot,
            'potentials_on_grid': self.plot_potentials_on_grid,
            'couplings_on_grid': self.plot_couplings_on_grid,
        }

    def plot(self, **kwargs):

        keys = self._define_main_keywords()
        what_keys = self._define_what_keywords()

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

    def full_residuals_plot(self, mlevels, path=None, fformat='png',
                            show=False, props=None):

        path = path or self.resid_path
        data = mlevels.out_data

        gen_color = self._get_color()
        gen_marker = self._get_marker()

        fillstyle = ['none', 'none', 'full', 'full']
        random.shuffle(fillstyle)

        props_def = {
            'xlabel': 'Energy',
            'ylabel': 'E_obs - E_calc',
            'title': 'Residuals',
        }

        # will change props values
        props = self._combine_props(props, props_def)

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
                'markersize': (6, 6),
                'markeredgewidth': (0.4, 0.4),
                'linewidth': (0.0, 0.0),
                'fillstyle': (fillstyle[0], fillstyle[1])
            }

            _, ax = plt.subplots()

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

            fig_path = os.path.join(path, 'residuals' + f'.{fformat}')
            plt.savefig(fig_path, format=fformat, dpi=250)  # transparent=True
            print(f'Created figure: {fig_path}', sep='')
            if show:
                plt.show()

    def plot_couplings_on_grid(self, mlevels, path=None, fformat='png',
                               show=False, props=None):

        path = path or self.func_path
        print(mlevels.fgrid)

        props_def = {
            'xlabel': 'Internuclear Distance',
            'ylabel': 'Interaction',
            'title': 'Couplings',
        }
        # will change props values
        props = self._combine_props(props, props_def)

        _, ax = plt.subplots()

        gen_color = self._get_color()

        for i, cpl in enumerate(mlevels.couplings):
            gridr = mlevels.rgrid * Const.bohr

            # will not work
            gridf_range = mlevels.fgrid[i*mlevels.ngrid:(i+1)*mlevels.ngrid]
            gridf = gridf_range * Const.hartree

            ax.plot(
                gridr,
                gridf,
                color=next(gen_color),
                markersize=0,
                linewidth=1.5,
            )

        ax.set_title(props['title'])
        ax.set_xlabel(props['xlabel'])
        ax.set_ylabel(props['ylabel'])
        ax.legend(mlevels.couplings, loc=0)

        plt.savefig(os.path.join(path, f'couplings.{fformat}'), format=fformat)

        if show:
            plt.show()

        # coupling_names = coupling_names or [ '' for _ in cols ]
        # for col, name in zip(cols, coupling_names):
        #     plt.plot(x, data[:,col])

    def plot_residuals_level(self, mlevels, path=None, fformat='png',
                             show=False, props=None):

        os.makedirs(self.res_path, exist_ok=True)
        path = path or self.res_path

        props_def = {
            'xlabel': 'J',
            'ylabel': 'E_obs - E_calc, cm-1',
            'title': 'v={}, state={}, isotop={}',
        }

        # will change props values
        props = self._combine_props(props, props_def)

        # if mlevels is not None:
        data = mlevels.out_data
        # else:
        #     data = np.loadtxt(evalues_file)

        # split data by v, then split by state, then extract
        # marker ranges and get the data for each isotope
        vsplit = np.split(data, np.where(np.diff(data[:, 1]))[0]+1)

        for v in range(0, len(vsplit)):
            ssplit = np.split(
                vsplit[v], np.where(np.diff(vsplit[v][:, -1]))[0]+1
            )

            for st in range(0, len(ssplit)):
                mrange = self._map_marker(ssplit[st][:, 7])

                for i, m in enumerate(mrange):
                    splitted_data = ssplit[st][m, :]

                    fdata = splitted_data[splitted_data[:, 6] == 0]
                    edata = splitted_data[splitted_data[:, 6] == 1]

                    fig, ax = plt.subplots()

                    ax.plot(
                        edata[:, 2],
                        -edata[:, 10],
                        color='blue',
                        marker='o',
                        markersize=10,
                        markeredgewidth=0.9,
                        fillstyle='none',
                        linewidth=0,
                    )

                    ax.plot(
                        fdata[:, 2],
                        -fdata[:, 10],
                        color='red',
                        marker='v',
                        markersize=8,
                        markeredgewidth=0.7,
                        fillstyle='full',
                        linewidth=0,
                    )

                    ax.set_xlabel(props['xlabel'])
                    ax.set_ylabel(props['ylabel'])
                    ax.legend(['e', 'f'], loc=0)
                    ax.set_title(props['title'].format(v, st, i))

                    figname = f'residual_{v}_{st}_{i}.{fformat}'
                    figpath = os.path.join(path, figname)
                    plt.savefig(figpath, format=fformat, dpi=100)
                    fig.tight_layout()
                    # plt.gcf().clear()

                    print(f'Created figure: {figpath}', sep='')
                    if show:
                        plt.show()

    def plot_potentials_on_grid(self, mlevels, path=None, fformat='png',
                                show=False, props=None):

        path = path or self.plot_path

        # to add x lim and y lim
        props_def = {
            'xlabel': 'Internuclear Distance',
            'ylabel': 'Potential energy',
            'title': 'Potentials',
        }

        # will change props values
        props = self._combine_props(props, props_def)

        _, ax = plt.subplots()

        unique_pfiles = set()
        for channel in mlevels.channels:
            unique_pfiles.add(channel.filep)

        # need to keep the unique file names in ordered data structue
        unique_pfiles = list(unique_pfiles)

        gen_color = self._get_color()

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

    def plot_potentials_points(self, files, path=None, fformat='png',
                               show=False, ipoints=None, xlim=None, ylim=None):
        # cannot not be called as option from plot()
        # only for plotting pointwise potentials

        _, ax = plt.subplots()

        path = path or self.plot_path

        gen_color = self._get_color()
        gen_marker = self._get_marker()
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

    def hcolormesh(self, mlevels, rows=None, cols=None, path=None,
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
        path = path or self.plot_path

        figpath = os.path.join(path, f'hcolormesh.{fformat}')
        plt.savefig(figpath, format=fformat)

        print(f'Created figure: {figpath}', sep='')

        if show:
            plt.show()

    def plot_wavefunctions(self, wffile, props=None, path=None,
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

        path = path or self.wavefunc_path

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

    def _combine_props(self, props, props_def):

        if props is None:
            return props_def
        else:
            for key, value in props_def.items():
                try:
                    props[key]
                except KeyError:
                    props[key] = value
            return props

    def _map_marker(self, markers):

        ranges = []
        for i in range(0, 60, 10):  # up to 6 isotopes
            res = np.where(np.logical_and(markers >= i, markers <= i+10))[0]
            if res.shape[0] > 0:
                ranges.append(res)

        return ranges

    def _get_color(self):

        bcolors = [
            '#1C03FF', '#4132C7', '#130B5B', '#10286B',
            '#355098', '#380AA2', '#3061B5', '#1C8FC9',
            '#00ABFF', '#10C7EC', '#10ECEC', '#BA5FFF'
        ]
        random.shuffle(bcolors)

        rcolors = [
            '#FA0A0A', '#D44141', '#B53232', '#EF262D',
            '#D50D50', '#FC697C', '#FC69E3', '#A11589',
            '#A8639C', '#DF211A', '#F3592A', '#FF8F44'
        ]
        random.shuffle(rcolors)

        gcolors = [
            '#00FF2B', '#3AA44C', '#165F22', '#54FC70',
            '#3EE22C', '#15730B', '#4EA644', '#65EE16',
            '#12DB7D', '#EAF215', '#DFF64D', '#F9F924'
        ]
        random.shuffle(gcolors)

        # colors = [
        #     '#e50000', '#9a0200', '#fe420f', '#f9bc08', '#9cbb04',
        #     '#419c03', '#69d84f', '#089404', '#06470c', '#01386a',
        #     '#448ee4', '#0d758f', '#1d5dec', '#0203e2', '#2000b1',
        #     '#380282', '#9a0eea', '#7e1e9c', '#4e0550', '#c20078',
        #     '#23c48b', '#fe86a4', '#ff474c', '#ffdf22', '#00122e'
        # ]
        # random.shuffle(colors)
        # for color in colors:
        #     yield color

        colors = []
        for color in zip(gcolors, bcolors, rcolors):
            g, b, r = color
            colors.append(g)
            colors.append(b)
            colors.append(r)

        for color in colors:
            yield color

    def _get_marker(self):

        markers = [
            'o', 'v', '<', '>', 's', 'p', '*',
            '+', 'x', 'D', 'd', 'X', '^'
        ]

        random.shuffle(markers)

        for marker in markers:
            yield marker
