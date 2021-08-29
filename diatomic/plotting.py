from os import makedirs as _makedirs
from os.path import join as _join
from random import shuffle as _shuffle
from scipy.interpolate import CubicSpline as _CubicSpline
import numpy as np
from utils import Utils, C_hartree, C_bohr
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib as mpl
import scipy.stats as st

mpl.use('TkAgg')

try:
    import seaborn as sns
    sns.set_style("white")
except ModuleNotFoundError:
    pass

__all__ = ['Plotting']


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
            'residuals2': cls.full_residuals_plot2,
            'residuals_hist': cls.full_residuals_hist,
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
    def plot_residuals_hist(cls, ml, path=None, fformat='png', save=True,
                            show=False, xlabel=None, ylabel=None,
                            bins=30, fname=None, with_sns=False, size=None):

        x = ml.out_data[:, 9]

        fig, ax = plt.subplots()
        if with_sns:
            sns.histplot(x, kde=True)
        else:
            n, bins, patches = ax.hist(
                x, density=True, bins=bins, label='Data'
            )
            mn, mx = plt.xlim()
            ax.set_xlim(mn, mx)
            kde_xs = np.linspace(mn, mx, 301)
            kde = st.gaussian_kde(x)
            ax.plot(kde_xs, kde.pdf(kde_xs), label='PDF')

        xlabel = r'Energy' or xlabel
        ylabel = 'Count' or ylabel
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        size = size or (6, 5)
        fig.set_size_inches(size)

        if save:
            path = path or Utils.get_plot_dir('resid_path')
            _makedirs(path, exist_ok=True)
            fname = fname or 'histogram'
            fig_path = _join(path, fname + f'.{fformat}')
            fig.savefig(fig_path, format=fformat, dpi=400)
            print(f'Created figure: {fig_path}', sep='')

        if show:
            plt.show()

    @classmethod
    def plot_residuals_by_isotopes(cls, ml, path=None, fformat='png',
                                   fname=None, xlabel=None, ylabel=None,
                                   show=False, msize=5, save=True, size=None):

        data = ml.out_data
        colors = [
            '#CF3F3F', '#4163A8', '#3FA561', '#DC8A56', '#B65CE6',
            '#D9CA5C', '#D165A9', '#99D469', '#171974'
        ]

        fig, ax = plt.subplots()

        for iso in ml.nisotopes:
            data_iso = data[(np.divmod(data[:, 6], 10)[0] + 1) == iso]

            ax.plot(
                data_iso[:, 8],
                data_iso[:, 9],
                color=colors[iso],
                marker='o',
                markersize=msize,
                markeredgewidth=0.8,
                linewidth=0.0,
                fillstyle='full',
                alpha=0.8,
                label=int(iso)
            )

        avrg_uncert = np.sum(data[:, 10]) / data.shape[0]
        hcolor = '#8B8E93'
        plt.axhline(
            y=avrg_uncert, color=hcolor, linestyle='dashed', linewidth=0.9
        )
        plt.axhline(
            y=-1.0*avrg_uncert, color=hcolor, linestyle='dashed', linewidth=0.9
        )

        xlabel = r'Energy' or xlabel
        yld = r'$\mathrm{E}_{\mathrm{calc}}$ - $\mathrm{E}_{\mathrm{obs}}$'
        ylabel = yld or ylabel
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(loc=0)
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.grid(which='major', alpha=0.4, linestyle='dashdot')
        ax.grid(which='minor', alpha=0.2, linestyle='dotted')

        size = size or (13, 7)
        fig.set_size_inches(size)

        if save:
            path = path or Utils.get_plot_dir('resid_path')
            _makedirs(path, exist_ok=True)
            fname = fname or 'residuals_ilabels'
            fig_path = _join(path, fname + f'.{fformat}')
            fig.savefig(fig_path, format=fformat, dpi=400)
            print(f'Created figure: {fig_path}', sep='')

        if show:
            plt.show()

    @classmethod
    def full_residuals_plot2(cls, ml, path=None, fformat='png', show=False):

        data = ml.out_data
        states = np.unique(data[:, -1])
        isotopes = np.unique(data[:, 6])

        colors = [
            '#AC4A4A', '#559ACF', '#56AD73', '#DC8A56', '#B65CE6',
            '#D9CA5C', '#D165A9', '#99D469', '#171974'
        ]

        for state in states:
            fig, ax = plt.subplots()
            count = 0
            for iso in isotopes:
                data_state = data[data[:, -1] == state]
                data_iso = data_state[data_state[:, 6] == iso]

                ax.plot(
                    data_iso[:, 8],
                    data_iso[:, 9],
                    color=colors[count],
                    marker='o',
                    markersize=6,
                    markeredgewidth=0.8,
                    linewidth=0.0,
                    fillstyle='full',
                    alpha=0.75,
                    label=int(iso)
                )
                count += 1

            avrg_uncert = np.sum(data[:, 10]) / data.shape[0]
            hcolor = '#8B8E93'
            plt.axhline(
                y=avrg_uncert, color=hcolor, linestyle='dashed', lw=0.9
            )
            plt.axhline(
                y=-1.0*avrg_uncert, color=hcolor, linestyle='dashed', lw=0.9
            )

            ax.set_xlabel('Energy')
            ax.set_ylabel('E_calc - E_obs')
            ax.legend(loc=0)

            ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
            ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
            ax.grid(which='major', alpha=0.4, linestyle='dashdot')
            ax.grid(which='minor', alpha=0.2, linestyle='dotted')

            path = path or Utils.get_plot_dir('resid_path')
            _makedirs(path, exist_ok=True)
            fig_title = 'residuals_slabels_' + str(int(state))
            fig.set_size_inches(13, 7)
            fig_path = _join(path, fig_title + f'.{fformat}')
            fig.savefig(fig_path, format=fformat, dpi=400)
            print(f'Created figure: {fig_path}', sep='')

            if show:
                plt.show()

    @classmethod
    def plot_residuals_by_states(cls, ml, path=None, fformat='png',
                                 fname=None, xlabel=None, ylabel=None,
                                 show=False, msize=5, save=True, size=None):

        data = ml.out_data
        states = np.unique(data[:, -1])

        colors = [
            '#AC4A4A', '#559ACF', '#56AD73', '#DC8A56', '#B65CE6',
            '#D9CA5C', '#D165A9', '#99D469', '#171974'
        ]

        fig, ax = plt.subplots()
        for state in states:
            data_state = data[data[:, -1] == state]

            ax.plot(
                data_state[:, 8],
                data_state[:, 9],
                color=colors[int(state)],
                marker='o',
                markersize=msize,
                markeredgewidth=0.8,
                linewidth=0.0,
                fillstyle='full',
                alpha=0.75,
                label=int(state)
            )

        avrg_uncert = np.sum(data[:, 10]) / data.shape[0]
        hcolor = '#8B8E93'
        plt.axhline(
            y=avrg_uncert, color=hcolor, linestyle='dashed', lw=0.9
        )
        plt.axhline(
            y=-1.0*avrg_uncert, color=hcolor, linestyle='dashed', lw=0.9
        )

        xlabel = r'Energy' or xlabel
        yld = r'$\mathrm{E}_{\mathrm{calc}}$ - $\mathrm{E}_{\mathrm{obs}}$'
        ylabel = yld or ylabel
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(loc=0)
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.grid(which='major', alpha=0.4, linestyle='dashdot')
        ax.grid(which='minor', alpha=0.2, linestyle='dotted')

        size = size or (13, 7)
        fig.set_size_inches(size)

        if save:
            path = path or Utils.get_plot_dir('resid_path')
            _makedirs(path, exist_ok=True)
            fname = fname or 'residuals_slabels_'
            fig_title = fname + str(int(state))
            fig_path = _join(path, fig_title + f'.{fformat}')
            fig.savefig(fig_path, format=fformat, dpi=400)
            print(f'Created figure: {fig_path}', sep='')

        if show:
            plt.show()

    @classmethod
    def plot_residuals_by_parity_labels(cls, ml, path=None, fformat='png',
                                        fname=None, xlabel=None, ylabel=None,
                                        show=False, msize=5, save=True, size=None):

        data = ml.out_data
        edata = data[data[:, 5] == 1]
        fdata = data[data[:, 5] == 0]

        line_styles = {
            'color': ('#B64A4A', '#4A73B5'),
            'marker': ('o', 'o'),
            'markersize': (msize, msize),
            'markeredgewidth': (0.85, 0.85),
            'linewidth': (0.0, 0.0),
            'fillstyle': ('full', 'full')
        }

        fig, ax = plt.subplots()
        ax.plot(
            edata[:, 8],
            edata[:, 9],
            color=line_styles['color'][0],
            marker=line_styles['marker'][0],
            markersize=line_styles['markersize'][0],
            markeredgewidth=line_styles['markeredgewidth'][0],
            linewidth=line_styles['linewidth'][0],
            fillstyle=line_styles['fillstyle'][0],
            alpha=0.85
        )

        ax.plot(
            fdata[:, 8],
            fdata[:, 9],
            color=line_styles['color'][1],
            marker=line_styles['marker'][1],
            markersize=line_styles['markersize'][1],
            markeredgewidth=line_styles['markeredgewidth'][1],
            linewidth=line_styles['linewidth'][1],
            fillstyle=line_styles['fillstyle'][1],
            alpha=0.85
        )

        # annotate average uncertanty
        avrg_uncert = np.sum(data[:, 10]) / data.shape[0]
        hcolor = '#8B8E93'
        plt.axhline(
            y=avrg_uncert, color=hcolor, linestyle='dashed', lw=0.9
        )
        plt.axhline(
            y=-avrg_uncert, color=hcolor, linestyle='dashed', lw=0.9
        )

        xlabel = r'Energy' or xlabel
        yld = r'$\mathrm{E}_{\mathrm{calc}}$ - $\mathrm{E}_{\mathrm{obs}}$'
        ylabel = yld or ylabel
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(['e-parity levels', 'f-parity levels'], loc=0)
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.grid(which='major', alpha=0.4, linestyle='dashdot')
        ax.grid(which='minor', alpha=0.2, linestyle='dotted')

        size = size or (13, 7)
        fig.set_size_inches(size)

        if save:
            path = path or Utils.get_plot_dir('resid_path')
            _makedirs(path, exist_ok=True)
            fname = fname or 'residuals_plabels'
            fig_path = _join(path, fname + f'.{fformat}')
            fig.savefig(fig_path, format=fformat, dpi=400)
            print(f'Created figure: {fig_path}', sep='')

        if show:
            plt.show()

    @classmethod
    def plot_couplings_on_grid(cls, mlevels, which=None, path=None,
                               show=False, save=False, fformat='png'):

        path = path or Utils.get_plot_dir('func_path')

        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        fig4, ax4 = plt.subplots()

        gen_color = cls._get_color()
        gen_marker = cls._get_marker()

        fgrid_cols = np.hstack((
            mlevels.rgrid[:, np.newaxis] * C_bohr,
            mlevels.fgrid.reshape(mlevels.ncp, mlevels.ngrid).T
        ))

        n = fgrid_cols.shape[0]

        leg1, leg2, leg3, leg4 = [], [], [], []

        for col in range(1, fgrid_cols.shape[1]):
            # TODO: chanage to real markers
            markers_on = [0, int(n/4), int(n/2), int(n/1.6), int(n/1.2), -1]

            if 'spin-orbit' in mlevels.couplings[col-1].coupling:
                ax1.plot(
                    fgrid_cols[:, 0],
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )
                ax1.set_xlabel('Internuclear Distance')
                ax1.set_ylabel('SO Interaction')
                # ax1.set_ylim(-800, -200)  # for nih
                # if ylim_so is not None:
                #     ax1.set_ylim(ylim_so[0], ylim_so[1])

                ax1.xaxis.set_minor_locator(tck.AutoMinorLocator())
                ax1.yaxis.set_minor_locator(tck.AutoMinorLocator())
                ax1.grid(which='major', alpha=0.5, linestyle='dashdot')
                ax1.grid(which='minor', alpha=0.25, linestyle='dotted')

                fig1.set_size_inches(8, 6)
                leg1.append(mlevels.couplings[col-1].interact)

            if 'LJ' in mlevels.couplings[col-1].coupling:
                ax2.plot(
                    fgrid_cols[:, 0],
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )
                ax2.set_xlabel('Internuclear Distance')
                ax2.set_ylabel('LJ Interaction')
                # ax2.set_ylim(0, 2) # for nih
                # if ylim_so is not None:
                #     ax2.set_ylim(ylim_so[0], ylim_so[1])

                ax2.xaxis.set_minor_locator(tck.AutoMinorLocator())
                ax2.yaxis.set_minor_locator(tck.AutoMinorLocator())
                ax2.grid(which='major', alpha=0.5, linestyle='dashdot')
                ax2.grid(which='minor', alpha=0.25, linestyle='dotted')

                fig2.set_size_inches(8, 6)
                leg2.append(mlevels.couplings[col-1].interact)

            if 'SJ' in mlevels.couplings[col-1].coupling:
                ax3.plot(
                    fgrid_cols[:, 0],
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )
                ax3.set_xlabel('Internuclear Distance')
                ax3.set_ylabel('SJ Interaction')
                # ax3.set_ylim(0, 3.5) # for nih
                # if ylim_so is not None:
                #     ax3.set_ylim(ylim_so[0], ylim_so[1])

                ax3.xaxis.set_minor_locator(tck.AutoMinorLocator())
                ax3.yaxis.set_minor_locator(tck.AutoMinorLocator())
                ax3.grid(which='major', alpha=0.5, linestyle='dashdot')
                ax3.grid(which='minor', alpha=0.25, linestyle='dotted')

                fig3.set_size_inches(8, 6)
                leg3.append(mlevels.couplings[col-1].interact)

            if 'LambdaD' in mlevels.couplings[col-1].coupling:
                ax4.plot(
                    fgrid_cols[:, 0],
                    fgrid_cols[:, col],
                    color=next(gen_color),
                    marker=next(gen_marker),
                    markersize=6.5,
                    markeredgewidth=0.7,
                    fillstyle='none',
                    linewidth=1.5,
                    markevery=markers_on
                )
                ax4.set_xlabel('Internuclear Distance')
                ax4.set_ylabel('2nd order Correction')
                # ax4.set_ylim(-0.1, 0.1) # for nih
                # if ylim_so is not None:
                #     ax4.set_ylim(ylim_so[0], ylim_so[1])

                ax4.xaxis.set_minor_locator(tck.AutoMinorLocator())
                ax4.yaxis.set_minor_locator(tck.AutoMinorLocator())
                ax4.grid(which='major', alpha=0.5, linestyle='dashdot')
                ax4.grid(which='minor', alpha=0.25, linestyle='dotted')

                fig4.set_size_inches(8, 6)
                leg4.append(mlevels.couplings[col-1].interact)

        ax1.legend(leg1, loc=0)
        ax2.legend(leg2, loc=0)
        ax3.legend(leg3, loc=0)
        ax4.legend(leg4, loc=0)

        if save:
            _makedirs(path, exist_ok=True)

            fig1_path = _join(path, f'SO.{fformat}')
            fig2_path = _join(path, f'LJ.{fformat}')
            fig3_path = _join(path, f'SJ.{fformat}')
            fig4_path = _join(path, f'LD.{fformat}')

            cls.save_figure(fig1, fig1_path, fformat, "tight")
            cls.save_figure(fig2, fig2_path, fformat, "tight")
            cls.save_figure(fig3, fig3_path, fformat, "tight")
            cls.save_figure(fig4, fig4_path, fformat, "tight")

        if show:
            plt.show()

    @classmethod
    def save_figure(self, fig, path, fformat, bbox):

        fig.savefig(path, format=fformat, bbox_inches=bbox, dpi=400)

    @classmethod
    def plot_residuals_level(cls, mlevels, path=None, fformat='png',
                             show=False, xlabel=None, ylabel=None,
                             save=True, size=None):

        res_path = Utils.get_plot_dir('res_path')
        _makedirs(res_path, exist_ok=True)
        path = path or res_path

        plt.rcParams.update({'font.size': 7})
        plt.rcParams.update({'figure.max_open_warning': 0})
        xlabel = xlabel or 'J'
        ylabel = ylabel or 'E_obs - E_calc [cm-1]'

        data = mlevels.out_data
        states = np.unique(data[:, -1])
        vibns = np.unique(data[:, 0])
        size = size or (6, 5)

        for v in vibns:
            datav = data[data[:, 0] == v]
            for state in states:
                fig, ax_data = plt.subplots()

                data_state = datav[datav[:, -1] == state]
                fdata = data_state[data_state[:, 5] == 0]
                edata = data_state[data_state[:, 5] == 1]

                ax_data.plot(
                    edata[:, 1], -edata[:, 9], color='#355896', marker='o',
                    markersize=4.0, markeredgecolor='#17376E', lw=0, alpha=0.8
                )

                ax_data.plot(
                    fdata[:, 1], -fdata[:, 9], color='#C43D3D', marker='X',
                    markersize=3.5, markeredgecolor='#903C3C', lw=0, alpha=0.8
                )

                ax_data.set_xlabel(xlabel)
                ax_data.set_ylabel(ylabel)
                ax_data.legend(
                    ['e', 'f'], loc=0, prop={'size': 6},
                    bbox_to_anchor=(1.1, 1.05)
                )
                ax_data.set_title(f'v={int(v)} | state={int(state)}', size=7)
                ax_data.xaxis.set_minor_locator(tck.AutoMinorLocator())
                ax_data.grid(which='major', alpha=1, linewidth=0.6)
                ax_data.grid(which='minor', alpha=1, linewidth=0.0)

                fig.set_size_inches(6, 5)

                if save:
                    figname = f'residual_{int(v)}_{int(state)}.{fformat}'
                    figpath = _join(path, figname)
                    fig.set_size_inches(size)
                    fig.savefig(figpath, format=fformat, dpi=400)
                    fig.tight_layout()

                    print(f'Created figure: {figpath}', sep='')

                if show:
                    plt.show()

    @classmethod
    def plot_residuals_level_subplots(cls, mlevels, path=None, fformat='png',
                                      show=False, save=True, size=None,
                                      nrows=1, ncols=2, isotope=None):

        data = mlevels.out_data
        states = np.unique(data[:, -1])
        vibns = np.unique(data[:, 0])

        fig, ax = plt.subplots(nrows, ncols, sharex=True, squeeze=False)

        row, col, count = 0, 0, 0

        for v in vibns:
            datav = data[data[:, 0] == v]
            for state in states:

                if col > ncols - 1:
                    col = 0
                    row += 1

                data_state = datav[datav[:, -1] == state]
                fdata = data_state[data_state[:, 5] == 0]
                edata = data_state[data_state[:, 5] == 1]

                if fdata.shape[0] == 0 or edata.shape[0] == 0:
                    continue

                ax[row, col].plot(
                    edata[:, 1], -edata[:, 9], color='#355896', marker='o',
                    markersize=5, lw=0, fillstyle='none', markeredgewidth=1.25
                )
                ax[row, col].plot(
                    fdata[:, 1], -fdata[:, 9], color='#C43D3D', marker='o',
                    markersize=5, lw=0, fillstyle='none', markeredgewidth=1.25
                )
                ax[row, col].grid(which='major', alpha=0.4)
                ax[row, col].grid(which='minor', alpha=0.2)
                ax[row, col].set_title(
                    f'v={int(v)} | state={int(state)}', size=10
                )
                ax[row, col].tick_params(
                    axis='both', which='major', labelsize=9
                )
                col += 1
                count += 1

        for i in range((ncols*nrows)-count):
            fig.delaxes(ax.flatten()[count+i])

        plt.rcParams.update({'font.size': 10})
        plt.rcParams.update({'figure.max_open_warning': 0})
        fig.legend(['e-level', 'f-level'])
        print(size)
        size = size or (14, 12)
        fig.set_size_inches(size)

        res_path = Utils.get_plot_dir('res_path')
        _makedirs(res_path, exist_ok=True)
        path = path or res_path

        if save:
            figname = f'residuals_1subplot.{fformat}'
            figpath = _join(path, figname)
            fig.savefig(figpath, format=fformat, dpi=400)
            fig.tight_layout()
            print(f'Created figure: {figpath}', sep='')

        if show:
            plt.show()

    @classmethod
    def plot_potentials_on_grid(cls, mlevels, path=None, fformat='png',
                                xlabel=None, ylabel=None, xlim=None, ylim=None,
                                show=False, save=True, size=None):

        unique_pfiles = set()
        for channel in mlevels.channels:
            unique_pfiles.add(channel.filep)
        unique_pfiles = list(unique_pfiles)

        gen_color = cls._get_color()
        fig, ax = plt.subplots()

        for i in range(mlevels.nch):
            gridr = mlevels.rgrid * C_bohr
            gridu_range = mlevels.ugrid[i*mlevels.ngrid:(i+1)*mlevels.ngrid]
            gridu = gridu_range * C_hartree
            ax.plot(gridr, gridu, color=next(gen_color), markersize=0, lw=1)

        xlim = xlim or ax.get_xlim()
        ylim = ylim or ax.get_ylim()
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        xlabel = xlabel or 'Internuclear Distance'
        ylabel = ylabel or 'Potential energy'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(unique_pfiles, loc=0)

        size = size or (6, 5)
        fig.set_size_inches(size)

        if save:
            path = path or Utils.get_plot_dir('plot_path')
            figpath = _join(path, f'potentials.{fformat}')
            plt.savefig(figpath, format=fformat, dpi=400)
            print(f'Created figure: {figpath}', sep='')

        if show:
            plt.show()

    @classmethod
    def plot_potentials_points(cls, files, ninter=50, path=None, fformat='png',
                               xlabel=None, ylabel=None, xlim=None, ylim=None,
                               show=False, save=True, size=None):

        gen_color = cls._get_color()
        gen_marker = cls._get_marker()
        size = size or (7, 5)
        _, ax = plt.subplots(figsize=size)

        for pfile in files:
            pdata = np.loadtxt(pfile, skiprows=1)
            x, y = pdata[:, 0], pdata[:, 1]
            x_interp = np.linspace(x[0], x[-1], ninter, endpoint=True)

            cs = _CubicSpline(x, y, bc_type='natural')
            y_interp = cs(x_interp)

            ax.plot(
                x_interp,
                y_interp,
                color=next(gen_color),
                marker=next(gen_marker),
                markersize=5.5,
                markeredgewidth=0.4,
                linewidth=1.2,
                fillstyle='none'
            )

        ax.set_title('Potentials')
        ax.set_xlabel('Internuclear distance')
        ax.set_ylabel('Potential energy')
        xlim = xlim or ax.get_xlim()
        ylim = ylim or ax.get_ylim()
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.legend(files, loc=0)

        if save:
            path = path or Utils.get_plot_dir('plot_path')
            figpath = _join(path, f'potential_points.{fformat}')
            plt.savefig(figpath, format=fformat, dpi=400)

            print(f'Created figure: {figpath}', sep='')

        if show:
            plt.show()

    @staticmethod
    def hcolormesh(mlevels, rows=None, cols=None, path=None,
                   fformat='png', show=False, save=True, size=None):

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
        size = size or (9, 9)
        fig, ax = plt.subplots(figsize=size)

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

        if save:
            # fig.tight_layout()
            path = path or Utils.get_plot_dir('plot_path')
            figpath = _join(path, f'hcolormesh.{fformat}')
            plt.savefig(figpath, format=fformat, dpi=400)
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

        _makedirs(path, exist_ok=True)

        if not subplots:

            fig_name = _join(path, 'wavefunctions' + f'.{fformat}')

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend(leg)
            plt.plot(igrid, wavefunc)
            plt.savefig(fig_name, format=fformat, dpi=400)

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
        _shuffle(bcolors)

        rcolors = [
            '#FA0A0A', '#D44141', '#B53232', '#EF262D',
            '#D50D50', '#FC697C', '#FC69E3', '#A11589',
            '#A8639C', '#DF211A', '#F3592A', '#FF8F44',
            '#C71585', '#D2691E', '#DC143C', '#800000'
        ]
        _shuffle(rcolors)

        gcolors = [
            '#00FF2B', '#3AA44C', '#165F22', '#54FC70',
            '#3EE22C', '#15730B', '#4EA644', '#65EE16',
            '#13987D', '#32B30B', '#77D11D', '#D1BE15',
            '#556B2F', '#2E8B57', '#008B8B', '#2E8B57'
        ]
        _shuffle(gcolors)

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

        _shuffle(markers)

        for marker in markers:
            yield marker
