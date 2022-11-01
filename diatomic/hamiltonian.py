import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline as _CubicSpline

from .identify import identify_levels
from .utils import Utils as _utils
from .interpolator import CSpline
from .operator import Operator


__all__ = ['Hamiltonian', 'PotEnr']


class Hamiltonian(Operator):
    """Construct the Hamiltonian operator, solve the radial
    Shrodinger equation and build a list of energy levels.

    Args:
        Operator (object): The parent Operator class
    """

    def __init__(self, objs, eig_decomp='lapack', lapack_driver='evr',
                 arpack_options=('LM', 6, None), is_weighted=False):

        """Construct the Hamiltonian object.

        Args:
            objs (list): List of objects
            eig_decomp (str, optional): Which package to be used for eigen
                decomposition, LAPACK or ARPACK. Defaults to 'lapack'.
            lapack_driver (str, optional): Defines which LAPACK driver
                should be used when eig_decomp is set to "lapack". Valid
                options are: "ev", "evd", "evr", "evx". syevr() procedure
                is called by default. Defaults to 'evr'.
            arpack_options (tuple, optional): A list of optional parameters
                about the arpack diagonalization. Defaults to ('LM', 6, None).
            is_weighted : (bool, optional). whether to apply Watson's weighting
                procedure. Defaults to False.
        """

        Operator.__init__(self, objs)

        pars, jnums = set(), set()
        for state in self.states:
            for symm in state.symmetry:
                pars.add(symm)

            for jrn in state.jqnumbers:
                jnums.add(jrn)

        self.pars = np.array(list(pars))
        self.jnums = np.array(list(jnums))

        if self.refj is not None:
            self._arange_jnums()

        # mapping arrays # self.Fy, self.Gy = grid.Fy, grid.Gy

        self.eig_decomp = eig_decomp
        self.lapack_driver = lapack_driver

        # which eigenvalues to use in the arpack procedure
        self.arpack_which = arpack_options[0]

        # in the arpack procedure k=6 by default
        self.arpack_k = arpack_options[1]

        # if sigma is set shift-invert mode is requested
        self.arpack_sigma = arpack_options[2]

        # map the diagonalization procedures
        self.diagonilize = {
            'lapack': self._lapack_eig_decomposition,
            'arpack': self._arpack_eig_decomposition
        }

        self.is_weighted = is_weighted
        # maximum number of fit parameters
        nmax_params = 200
        # maximum number of computed levels
        self.nmax_levels = 80000
        # matrix with the spline S functions
        self.sk_grid = np.zeros((self.msize, nmax_params))

        self.kin_enr = KinEnr(objs)

    def _arange_jnums(self):

        self.jnums = np.delete(self.jnums, np.where(self.jnums == self.refj))
        self.jnums = np.insert(self.jnums, 0, self.refj)

    def get_energy_levels(self, ident_option=1):
        """Generate a list of identified energy levels

        Args:
            ident_option (int, optional): What method for identification
            to use. Defaults to 1.

        Raises:
            AttributeError: In case of missing observations

        Returns:
            array: The generated list of levels
        """

        self.ident_option = ident_option if ident_option in [0, 1] else 1

        if self.exp_data is not None:
            # edatas = self.exp_data[self.exp_data[:, 7].argsort()]
            # cdatas = self.calc_data[self.calc_data[:, -4].argsort()]
            # cdatas = self.calc_data[self.calc_data[:, 1].argsort()]
            # cdatas = self.calc_data[self.calc_data[:, -3].argsort()]
            # cdatas = self.calc_data[self.calc_data[:, 2].argsort()]
            # cdatas = self.calc_data[self.calc_data[:, -4].argsort()]

            # edatas = self.exp_data[self.exp_data[:, 3].argsort()]
            # edatas = self.exp_data[self.exp_data[:, 1].argsort()]
            # edatas = self.exp_data[self.exp_data[:, 2].argsort()]
            # edatas = self.exp_data[self.exp_data[:, 7].argsort()]

            self.out_data = identify_levels(
                np.ascontiguousarray(self.calc_data),
                np.ascontiguousarray(self.exp_data),
                self.nch,
                ident_option)
        else:
            raise AttributeError(
                'Experimental energy terms should be provided.')

        return self.out_data

    def solve(self, ops, ypar=None, energy_range_index=None, energy_range_value=None):
        """A general function which solves the coupled system of Shcrodinger equations

        Args:
            energy_range_index (tuple, optional): Specifies a subset of
                energies by providing two indices: the initial and the final
                index of the coumputed eigenvalues. Defaults to None.
            energy_range_value (tuple, optional): Specifies a subset of
                energies by providing two energy values: the start and the end
                value of the coumputed eigenvalues. Defaults to None.

        Returns:
            tuple: A tuple of the computed eigenvalues and eigenvectors
        """
        self.ops = ops
        self.energy_range_value = energy_range_value
        self.energy_range_index = energy_range_index

        if ypar is None:
            for operator in self.ops:
                operator.calculate_radial_function_on_grid(self.params)
        else:
            for operator in self.ops:
                operator.calculate_radial_function_on_grid(ypar)

        # count the total number of eigenvalues
        self.ecount = 0

        evalues_all = np.zeros((self.nmax_levels, 9+self.nch))
        evecs_all = np.zeros((self.msize, self.nmax_levels))

        for _, iso in enumerate(self.niso):
            tmatrix = self.kin_enr.calculate_kinetic_energy(iso)

            for par in self.pars:

                shift_enr = [0.0]  # pass by reference; make it self?
                for jrotn in self.jnums:

                    # build and diagonalize the hamiltonian matrix
                    hmatrix = self.build_hamiltonian(iso, par, jrotn, tmatrix)
                    evalues, evecs = self.diagonilize[self.eig_decomp](hmatrix)
                    nevalues = evalues.shape[0]

                    # arange and store the quantum numbers & labels for levels
                    levels = self._arange_levels(
                        jrotn, par, iso, evalues, evecs, shift_enr=shift_enr)
                    evalues_all[self.ecount:nevalues+self.ecount, :] = levels

                    evecs_all[:, self.ecount:nevalues+self.ecount] = evecs
                    self.ecount += nevalues

        evalues_all = evalues_all[:self.ecount, :]
        evecs_all = evecs_all[:, :self.ecount]

        self.hamiltonian_matrix = hmatrix

        self.calc_data, self.evecs_matrix = evalues_all, evecs_all

        return evalues_all, evecs_all

    def wrapper_calculate_diatomic_levels(self, ypar):
        """A wrapper procedure for calculating the energy levels
        which is called through the fitting routine

        Args:
            ypar (array): The fitting parameters

        Returns:
            array: The computed energy levels
        """
        # self.interpolate_functions(ypar)
        self.solve(self.ops, ypar=ypar, energy_range_index=self.energy_range_index,
                   energy_range_value=self.energy_range_value)
        return self.get_energy_levels(self.ident_option)

    def extract_terms_in_range(self, uenergy=None, lenergy=None, usymm=None,
                               lsymm=None, uj=None, lj=None, ustate=None,
                               lstate=None):

        # TODO: it may happen that the number of levels within the specified
        # ranges will vary during the fit

        eind, jind, pind, sind = 1, 2, 3, -4

        # filter by energy
        self.uenergy = uenergy or self.calc_data[0, eind]
        self.lenergy = lenergy or self.calc_data[0, eind]
        calc_uterms = self.calc_data[
            (self.calc_data[:, eind] >= self.uenergy[0]) &
            (self.calc_data[:, eind] <= self.uenergy[1])]
        calc_lterms = self.calc_data[
            (self.calc_data[:, eind] >= self.lenergy[0]) &
            (self.calc_data[:, eind] <= self.lenergy[1])]

        # filter by J
        self.uj = uj or (np.amin(calc_uterms[:, jind]),
                         np.amax(calc_uterms[:, jind]))
        self.lj = lj or (np.amin(calc_lterms[:, jind]),
                         np.amax(calc_lterms[:, jind]))
        calc_uterms = calc_uterms[
            (calc_uterms[:, jind] >= self.uj[0]) &
            (calc_uterms[:, jind] <= self.uj[1])]
        calc_lterms = calc_lterms[
            (calc_lterms[:, jind] >= self.lj[0]) &
            (calc_lterms[:, jind] <= self.lj[1])]

        # filter by symmetry
        self.usymm = usymm or (0, 1)
        self.lsymm = lsymm or (0, 1)
        calc_uterms = calc_uterms[
            (np.in1d(calc_uterms[:, pind], self.usymm[0])) |
            (np.in1d(calc_uterms[:, pind], self.usymm[1]))]
        calc_lterms = calc_lterms[
            (np.in1d(calc_lterms[:, pind], self.lsymm[0])) |
            (np.in1d(calc_lterms[:, pind], self.lsymm[1]))]

        # filter by state - use only for the case of single state
        self.ustate = ustate
        self.lstate = lstate
        if self.ustate is not None:
            calc_uterms = calc_uterms[calc_uterms[:, sind] == ustate]
        if self.lstate is not None:
            calc_lterms = calc_lterms[calc_lterms[:, sind] == lstate]

        return calc_uterms, calc_lterms

    def build_hamiltonian(self, iso, par, jrotn, tmatrix):
        """Builds the Hamiltonian matrix as a sum of three matricies:
        the kinetic energy, the potential energy and the couplings

        Args:
            iso (int): isotopologue number
            par (int): e- or f-symmetry
            jrotn (float): the rotational quantum number
            tm (array): the kinetic energy matrix

        Returns:
            array: the final constructed Hamiltonian matrix
        """

        hmatrix = np.zeros((self.msize, self.msize))
        hmatrix += tmatrix

        for opr in self.ops:
            ch1, ch2 = opr.istate1, opr.istate2
            ri1, ri2 = (ch1-1)*self.ngrid, ch1*self.ngrid  # row indices
            ci1, ci2 = (ch2-1)*self.ngrid, ch2*self.ngrid  # col indices

            if par in opr.state1.symmetry or par in opr.state2.symmetry:
                matrix_block = opr.calculate_matrix_elements(jrotn, par, iso)
                hmatrix[ri1:ri2, ci1:ci2] += matrix_block
            else:
                hmatrix[ri1:ri2, ci1:ci2] = np.zeros((self.ngrid, self.ngrid))

        return hmatrix

    def _lapack_eig_decomposition(self, hmatrix):
        """Diagonilizes the Hamiltonian matrix through the scipy
        eigh() procedure from the LAPACK package

        Args:
            hmatrix (array): The Hamiltonian

        Returns:
            tuple: A tuple of eigenvalues and eigenvectors
        """

        subset_value, subset_index = None, None

        if self.energy_range_index:
            emini = self.energy_range_index[0]
            emaxi = self.energy_range_index[-1]
            subset_index = (emini, emaxi-1)

        elif self.energy_range_value:
            eminv = self.energy_range_value[0] / _utils.C_hartree
            emaxv = self.energy_range_value[-1] / _utils.C_hartree
            subset_value = (eminv, emaxv)
        else:
            subset_value = (-np.inf, np.inf)

        evalues, evecs = sp.linalg.eigh(
            hmatrix,
            eigvals_only=False,
            overwrite_a=True,
            lower=False,
            subset_by_index=subset_index,
            subset_by_value=subset_value,
            driver=self.lapack_driver.lower(),
            check_finite=False
        )

        return evalues, evecs

    def eigsh(self):
        raise NotImplementedError()

    def eigh(self):
        raise NotImplementedError()

    def _arpack_eig_decomposition(self, hmatrix):
        """Diagonilizes the Hamiltonian matrix through the scipy sparse eigsh()
        procedure from the ARPACK package

        Args:
            hmatrix (array): The Hamiltonian matrix

        Returns:
            tuple: A tuple of eigenvalues and eigenvectors

        Notes:
            ARAPCK procedure is the most efficient and approppriate for finding
            the largest eigenvalues of a sparse matrix. If the smallest
            eigenvalues are required then it is recommended to use a
            shift-invert mode. It transforms the eigenvalue problem to
            an eqivalent problem with shifted eigenvalues in which the
            small eigenvalues u become the large eigenvalues v: v = 1 / u
        """

        evalues, evecs = sp.sparse.linalg.eigsh(
            hmatrix,
            k=self.arpack_k,
            which=self.arpack_which.upper(),
            sigma=self.arpack_sigma,
            return_eigenvectors=True)

        return evalues, evecs

    def _arange_levels(self, jrotn, par, iso, evalues, evecs, shift_enr=[0.0]):
        """ Given a set of good quantum numbers, labels and the computed
        eigenvalues, this function will compute additional information
        about the energy levels and arange the entire information for
        each level in one matrix

        Args:
            jrotn (float): The rotational quantum number
            par (int): e- or f-symmetry
            iso (int): isotopolgue number
            evalues (array): The array of computed eigenvalues
            evecs (array): The computed eigenvectors
            shift_enr (list, optional): The shift energy. Defaults to [0.0].

        Returns:
            array: The energy levels with the corresponding information
                aranged by columns
        """

        ids = np.arange(self.ecount+1, evalues.shape[0]+self.ecount+1)

        if self.refj and jrotn == self.jnums[0]:
            shift_enr[0] = evalues[0]
        elif self.ref_enr:
            shift_enr[0] = self.ref_enr

        evalues_shifted = (evalues - shift_enr[0]) * _utils.C_hartree

        ccoefs = self._get_coupling_coefficients(evecs)
        states = np.argmax(ccoefs, axis=1) + 1
        vibnums = self._assign_vibrational_number(states)
        lambdas = self._get_lambda_values(states)
        omegas = self._get_omega_values(states)

        levels = np.column_stack((
            ids, evalues_shifted,
            np.full((evalues_shifted.shape[0],), jrotn),
            np.full((evalues_shifted.shape[0],), par),
            np.full((evalues_shifted.shape[0],), iso),
            ccoefs, states, vibnums, lambdas, omegas))

        # shift energies
        # for i, ch in enumerate(self.states):
        #     if ch.shift_by != 0.0:
        #         inds = np.where([levels[:, -4] == i+1])[1]
        #         levels[inds, 1] += ch.potential.shift_by * (-1)

        return levels

    def _get_coupling_coefficients(self, evecs):
        """Compute the coupling coeffcients which can be further use
        for assigning a state based on the largest coefficient

        Args:
            evecs (array): The eigenvectors array

        Returns:
            array: The computed coeffcients
        """

        ccoefs = np.zeros((evecs.shape[1], self.nch))
        evecs = np.power(evecs, 2)

        for k in range(0, self.nch):
            ccoefs[:, k] = evecs[k*self.ngrid:(k+1)*self.ngrid, :].sum(axis=0)

        return ccoefs

    def _assign_vibrational_number(self, states):

        vibns = np.empty(states.shape[0])
        vibns_count = np.empty(self.nch)
        vibns_count.fill(-1)

        for i, state in enumerate(states):
            vibns_count[state-1] += 1
            vibns[i] = vibns_count[state-1]

        return vibns

    def _get_lambda_values(self, states):

        return np.fromiter(
            map(lambda x: self.states[x-1]._lambda, np.int64(states)),
            dtype=np.int64)

    def _get_omega_values(self, states):

        return np.fromiter(
            map(lambda x: self.states[x-1].omega, np.int64(states)),
            dtype=np.float64)

    def add_second_order_correction(self, state, qe=0., qf=0.):

        finds = np.where(
            (self.calc_data[:, 3] == 0) & (self.calc_data[:, -4] == state))
        einds = np.where(
            (self.calc_data[:, 3] == 1) & (self.calc_data[:, -4] == state))

        jjf = self.calc_data[finds, 2] * (self.calc_data[finds, 2] + 1.0)
        self.calc_data[finds, 1] = self.calc_data[finds, 1] + qf*jjf

        jje = self.calc_data[einds, 2] * (self.calc_data[einds, 2] + 1.0)
        self.calc_data[einds, 1] = self.calc_data[einds, 1] + qe*jje

    def save_energy_list(self, levels=None, stats=None, filename=None):

        out_levels = self.out_data if levels is None else levels

        if out_levels.shape[0] == 0:
            raise ValueError('Empty list of levels!')

        if stats is None:
            stats = _utils.calculate_stats(out_levels[:, 8], out_levels[:, 7],
                                           out_levels[:, 10], self.is_weighted)

        # add id column
        out_levels = np.c_[np.arange(1, out_levels.shape[0]+1), out_levels]

        coef = 'coef'
        coef_labels = ''.join(
            map(lambda x: f'{coef+str(x):^8}', range(1, self.nch+1)))

        labels = ('n', 'v', 'J', 'omega', 'sigma', 'lambda',
                  'symmetry', 'marker', 'Ecalc', 'Eexp',
                  'delta', 'unc', coef_labels, 'state')

        header = (f'{labels[0]:^11}{labels[1]:<7}{labels[2]:<6}{labels[3]:<8}'
                  f'{labels[4]:<6}{labels[5]:<7}{labels[6]:<9}{labels[7]:<10}'
                  f'{labels[8]:<14}{labels[9]:<15}{labels[10]:<12}'
                  f'{labels[11]:<7}{labels[12]:>15}{labels[13]}')

        fmt = ['%7d', '%5d', '%7.1f', '%7.1f', '%7.1f', '%5d', '%6d',
               '%6d', '%15.6f', '%14.6f', '%12.6f', '%9.4f']
        fmt += self.nch*['%7.3f'] + ['%5d']

        footer = _utils.print_stats(stats)

        file_energies = f'{self.mol_name}_energy_levels.dat'
        out_fname = filename or file_energies

        np.savetxt(
            out_fname, out_levels, header=header, footer=footer, fmt=fmt)

    def save_full_energy_list(self, calc_energies=None, filename=None):
        """Stores the complete list of computed energy levels in external file

        Args:
            calc_energies (array): the array containing the energy levels
            filename (str, optional): The name of the file. Defaults to None.
        """

        cal_data = self.calc_data if calc_energies is None else calc_energies

        if cal_data.shape[0] == 0:
            raise ValueError('Empty list of eigenvalues!')

        nrows = cal_data.shape[1]
        cols = [0, -3] + list(range(1, nrows-3)) + [-2, -1]
        # calc_data_out = self.calc_data[:, cols]
        cal_data = cal_data[:, cols]

        coef = 'coef'
        coef_labels = ''.join(
            map(lambda x: f'{coef+str(x):^10}', range(1, self.nch+1)))

        labels = ('n', 'v', 'Ecalc', 'J', 'symmetry', ' marker',
                  coef_labels, 'state', 'lambda', 'omega')

        header = (f'{labels[0]:^10}{labels[1]:<9}{labels[2]:<15}{labels[3]:<4}'
                  f'{labels[4]:<7}{labels[5]:<7}{labels[6]}{labels[7]:<7}'
                  f'{labels[8]:<9}{labels[9]}')

        fmt = ['%7.1d', '%5.1d', '%16.6f', '%7.1f', '%5.1d', '%7.1d']
        fmt += self.nch*['%9.3f'] + 2*['%6.1d'] + ['%8.1f']

        file_energies = f'{self.mol_name}_eigenenergies.dat'
        output_file = filename or file_energies
        np.savetxt(output_file, cal_data, header=header, fmt=fmt)

    def _sort_full_energy_list(self, cols):

        out = [self.calc_data[:, col] for col in reversed(cols)]
        self.calc_data = self.calc_data[np.lexsort(tuple(out))]

    def sort_full_energy_list(self, cols=[]):

        self._sort_full_energy_list(cols)
        self.save_full_energy_list()

    def _sort_energy_levels(self, cols):

        cols = [i-2 for i in cols]
        out = [self.out_data[:, col] for col in reversed(cols)]
        self.out_data = self.out_data[np.lexsort(tuple(out))]

    def sort_energy_levels(self, cols=[]):

        self._sort_energy_levels(cols)
        self.save_energy_list()

    def get_predicted_data(self):
        """Get the complete list of computed energy levels

        Returns:
            array: the computed energy levels
        """

        return self.calc_data

    def get_output_data(self):
        """Get the full list of energy levels

        Returns:
            array: The computed list of energy levels
        """

        return np.c_[np.arange(1, self.out_data.shape[0]+1), self.out_data]

    def get_hamiltonian(self):
        """Get the Hamiltonian matrix for the last computed J, symmetry
        and isotopologue

        Returns:
            array: The Hamiltonian matrix
        """

        return self.hamiltonian_matrix

    def interp_wavefunc(self):

        ninter = 600
        rgrid, rstep = self.rgrid * _utils.C_bohr, self.rstep * _utils.C_bohr
        igrid, istep = np.linspace(self.rmin, self.rmax, num=ninter,
                                   endpoint=True, retstep=True)

        # sinc_matrix = np.zeros((self.ngrid, igrid.shape[0]))
        res = np.zeros(ninter)

        ncoefs = self.evecs_matrix[:2*self.ngrid, 0]

        k = 0
        for j in range(igrid.shape[0]):
            for i in range(2*self.ngrid):
                if k >= self.ngrid:
                    k = 0

                argi = (igrid[j] - rgrid[k]) / rstep
                sinc = np.sinc(argi)
                # sinc_matrix[i, :] = sinc
                res[j] += ncoefs[i] * sinc

        return igrid, res

    def plot_radial_functions(self, ops, show=False, fname=None, subplots=None,
                              fsize=None, xlim=None, ylim=None):
        """Plot the radial functions on grid

        Args:
            ops (list): List of operators whose radial functions to be plotted
            show (bool, optional): Thether to show the plot. Defaults to False.
            fname (str, optional): The file name of the plot. Defaults to None.
            subplots (tuple, optional): The number of subplots in x and y dim.
            Defaults to None.

        Raises:
            ValueError: If the subplots format is incorrect.
        """

        plt.style.use('seaborn-paper')
        plt.style.context('seaborn-paper')
        fsize = fsize or (6, 5)

        if not subplots:
            fig, ax = plt.subplots(1, 1, figsize=fsize, constrained_layout=True)

            for opr in ops:
                ax.plot(opr.rgrid, opr.ygrid, lw=1.2, label=opr.label)
                ax.tick_params(axis='both', direction="in", which='major', labelsize=11)
                ax.legend(loc='best', fontsize=11)
                if xlim is not None:
                    ax.set_xlim(xlim[0], xlim[1])
                if ylim is not None:
                    ax.set_ylim(ylim[0], ylim[1])
        else:
            try:
                nrows, ncols = subplots[0], subplots[1]
            except (ValueError, TypeError) as exc:
                raise ValueError('Subplots should be a tuple of 2 values.') from exc

            fig, axs = plt.subplots(nrows, ncols, figsize=fsize,
                                    constrained_layout=True)
            axs = axs.flatten()
            for ax, opr in zip(axs, ops):
                ax.plot(opr.rgrid, opr.ygrid, lw=1.2)
                ax.set_title(label=opr.label)
                ax.tick_params(axis='both', direction="in", which='major', labelsize=10)
                if xlim is not None:
                    ax.set_xlim(xlim[0], xlim[1])
                if ylim is not None:
                    ax.set_ylim(ylim[0], ylim[1])

        fig.supxlabel('R ' + '(Angstrom)', fontsize=11)
        fig.supylabel('Radial function', fontsize=11)

        if show:
            plt.show()

        if fname is not None:
            self._save_plot(fig, fname)

    def plot_hamiltonian_colormesh(self, rows=None, cols=None, show=False,
                                   fname=None, fsize=None):
        """Create a colormesh of the Hamiltonian matrix

        Args:
            rows (tuple, optional): The row indices. Defaults to None.
            cols (tuple, optional): The col indices. Defaults to None.
            show (bool, optional): Whether to show the image. Defaults to False.
            fname (str, optional): The file name of the created figure.
            Defaults to None.

        Remarks:
            1. The parameters rows and cols are defined by two integer numbers -
            the first and the last row and column, respectively of the submatrix
            whcih should be plotted.
            2. vmin and vmax set the normalization range.
            By default scale scalar data to [0, 1] range.
        """

        # TODO: The Hamiltonian is upper diagonal matrix - convert it to full
        hmat = np.copy(self.hamiltonian_matrix)
        hmat[hmat == 0.0] = np.nan

        rows = rows or (0, hmat.shape[0])
        cols = cols or (0, hmat.shape[1])
        fsize = fsize or (9, 7)
        fig, ax = plt.subplots(figsize=fsize, constrained_layout=True)

        # im = ax.pcolormesh(hmat[rows[0]:rows[1], cols[0]:cols[1]],
        # vmin=1.0e-8, vmax=0.01, cmap='RdBu_r')

        img = ax.matshow(hmat[rows[0]:rows[1], cols[0]:cols[1]],
                        vmin=1.0e-8, vmax=0.01, cmap='RdBu_r',
                        aspect='auto', interpolation=None)

        fig.colorbar(img)
        ax.set_title(
            f'Hamiltonian matrix, rows={rows[0]}:{rows[1]}, '
            f'cols={cols[0]}:{cols[1]}')

        if show:
            plt.show()

        if fname is not None:
            self._save_plot(fig, fname)

    def plot_eigenfunctions(self, nlevels, show=False, fname=None,
                            subplots=None, fsize=None):

        plt.style.use('seaborn-paper')
        plt.style.context('seaborn-paper')

        fsize = fsize or (6, 5)
        subplots = subplots or (1, 1)
        try:
            _ = iter(nlevels)
        except TypeError:
            nlevels = list(nlevels)

        fontsize = 12
        try:
            nrows, ncols = subplots[0], subplots[1]
        except (ValueError, TypeError) as exc:
            raise ValueError('Subplots should be a tuple of 2 values.') from exc

        fig, axs = plt.subplots(nrows, ncols, figsize=fsize, sharex=True,
                                constrained_layout=True)

        try:
            _ = iter(axs)
            for ax, nlevel in zip(axs, nlevels):
                evec = self.evecs_matrix[:, nlevel-1]  # substrct 1 since the counting starts from 0
                ax.plot(self.rgrid, evec, lw=1.2, label=str(nlevel))
                ax.set_title(label=f'level {str(nlevel)}')
                ax.tick_params(axis='both', direction="in", which='major', labelsize=10)
        except TypeError:
            for nlevel in nlevels:
                evec = self.evecs_matrix[:, nlevel-1]  # substrct 1 since the counting starts from 0
                axs.plot(self.rgrid, evec, lw=1.2, label=str(nlevel))
                axs.tick_params(axis='both', direction="in", which='major', labelsize=10)
                axs.legend(loc='best', fontsize=fontsize)

        fig.supxlabel('R ' + '(Angstrom)', fontsize=fontsize)
        fig.supylabel('Wavefunction', fontsize=fontsize)

        if show:
            plt.show()

        if fname is not None:
            self._save_plot(fig, fname)

    def plot_levels_hist(self, data=None, show=False, fname=None,
                         fsize=None, bins=30, sns=False):

        plt.style.use('seaborn-paper')
        plt.style.context('seaborn-paper')
        fsize = fsize or (6, 5)
        fig, ax = plt.subplots(1, 1, figsize=fsize, constrained_layout=True)

        if data is None:
            data = self.out_data[:, 9]

        if sns:
            try:
                import seaborn as sns
                sns.set_style("white")
                sns.histplot(data, kde=True)
            except ModuleNotFoundError:
                print('Seaborn package is not installed.')
        else:
            n, bins, patches = ax.hist(data, density=True, bins=bins)
            xmin, xmax = plt.xlim()
            ax.set_xlim(xmin, xmax)
            kde_xs = np.linspace(xmin, xmax, 301)
            kde = sp.stats.gaussian_kde(data)
            ax.plot(kde_xs, kde.pdf(kde_xs), label='PDF')

        fig.supxlabel('Energy', fontsize=12)
        fig.supylabel('Count', fontsize=12)

        if fname is not None:
            self._save_plot(fig, fname)

        if show:
            plt.show()

    def _save_plot(self, fig, fname):

        try:
            frmt = os.path.basename(fname).split('.')[1]
        except IndexError:
            frmt = 'pdf'
        fig.savefig(fname, format=frmt, dpi=300)


class KinEnr(Operator):
    """Calculate the kinetic energy operator
    """
    def __init__(self, objs):

        Operator.__init__(self, objs)

        self.T = None

    def calculate_kinetic_energy(self, iso):
        """Call a specific function for the calculation of the kinetic
        energy operator depending on the chosen solver method

        Args:
            iso (int): The isotopologue number

        Raises:
            ValueError: If the name of the chosen solver method is not
                present in the list of allowed values

        Returns:
            array: The computed kinetic energy matrix
        """

        mass = self.masses[iso-1]

        if self.solver == 'sinc':
            return self._calculate_kinetic_energy_sinc(mass)

        elif self.solver == 'fourier':
            return self._calculate_kinetic_energy_fourier(mass)

        elif self.solver == 'fd5':
            return self._calculate_kinetic_energy_fd5(mass)

        else:
            raise ValueError(
                f'Error: {self.solver} is not allowed solver method.')

    def _calculate_kinetic_energy_fourier(self, mass):
        """Calculate the kinetic energy operator using the Fourir basis

        Args:
            mass (double): The mass value

        Returns:
            array: The computed kinetic energy matrix
        """
        T = np.zeros((self.msize, self.msize))

        length = self.rmax - self.rmin
        n2 = (self.ngrid**2 + 2.0) / 6.0
        planck_const = 1.0  # in au
        ml = 4.0 * mass * length**2
        h2 = (2.0 * np.pi * planck_const)**2

        ij_grid = np.mgrid[:self.ngrid, :self.ngrid]
        ij_mtx = ij_grid[0] - ij_grid[1]  # TODO: check this
        diag_inds = np.diag_indices(self.ngrid)
        sinf_mtx = np.sin((ij_mtx * np.pi) / self.ngrid)

        # set diagonal to some temporary nonzero value
        sinf_mtx[diag_inds] = 1.0

        T[:self.ngrid, :self.ngrid] = \
            (np.power(-1.0, ij_mtx) * h2) / (ml * np.power(sinf_mtx, 2))

        T[diag_inds] = (h2 / ml) * n2

        for chn in range(2, self.nch+1):
            ind1, ind2 = (chn-1)*self.ngrid, chn*self.ngrid
            T[ind1:ind2, ind1:ind2] = T[:self.ngrid, :self.ngrid]

        self.T = T

        return T

    def _calculate_kinetic_energy_sinc(self, mass):
        """Calculate the kinetic energy operator using the sinc basis

        Args:
            mass (double): The mass value

        Returns:
            array: The computed kinetic energy matrix
        """
        T = np.zeros((self.msize, self.msize))

        planck_const = 1.0  # in au
        pi2 = np.pi ** 2
        hstep = (self.rmax - self.rmin) / (self.ngrid-1)
        h2 = (planck_const ** 2) / (2.0 * mass)

        ij_grid = np.mgrid[:self.ngrid, :self.ngrid]
        ij_diff = ij_grid[0] - ij_grid[1]
        ij_sum = ij_grid[0] + ij_grid[1]
        diag_inds = np.diag_indices(self.ngrid)

        # set diagonal to some temporary nonzero value
        ij_diff[diag_inds] = 1.0

        T[:self.ngrid, :self.ngrid] = \
            (2.0 * np.power(-1.0, ij_sum)) / np.power(ij_diff, 2)
        T[diag_inds] = pi2 / 3.0
        T = T * (h2 / hstep**2)

        for chn in range(2, self.nch+1):
            ind1, ind2 = (chn-1)*self.ngrid, chn*self.ngrid
            T[ind1:ind2, ind1:ind2] = T[:self.ngrid, :self.ngrid]

        # T[di] = T[di] - (h2 * self.Fy)

        self.T = T

        return T

    def _calculate_kinetic_energy_fd5(self, mass):
        """Calculate the kinetic energy operator using the finite difference method

        Args:
            mass (double): The mass value

        Returns:
            array: The computed kinetic energy matrix
        """

        T = np.zeros((self.msize, self.msize))

        # the first and last 2 eigenvalues are wrong
        planck_const = 1.0  # in au
        h2 = (planck_const ** 2) / (2.0 * mass)
        step = (self.rmax - self.rmin) / (self.ngrid-1)

        diag0 = np.empty(self.ngrid)
        diag0.fill(5.0/2.0)

        diag1 = np.empty(self.ngrid-1)
        diag1.fill(-4.0/3.0)

        diag2 = np.empty(self.ngrid-2)
        diag2.fill(1.0/12.0)

        T[:self.ngrid, :self.ngrid] = \
            (h2/(step**2)) * self._five_diagonal_matrix(diag0, diag1, diag2)

        corner_coef = 29.0 / 12.0
        T[0, 0] = corner_coef
        T[self.ngrid-1, self.ngrid-1] = corner_coef

        for chn in range(2, self.nch+1):
            ind1, ind2 = (chn-1)*self.ngrid, chn*self.ngrid
            T[ind1:ind2, ind1:ind2] = T[:self.ngrid, :self.ngrid]

        self.T = T

        return T

    def _five_diagonal_matrix(self, a, b, c, k1=0, k2=1, k3=2):

        return np.diag(a, k1) + np.diag(b, k2) + \
               np.diag(b, -k2) + np.diag(c, k3) + \
               np.diag(c, -k3)


class PotEnr(Operator):
    """Calculate the Potential energy operator

    Args:
        Operator (object): The parent class
    """

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)

        self.set_states()

        if self.istate1 != self.istate2:
            raise ValueError(
                "The states in 'pair_states' should coincide for 'PotEnr'.")

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        mass = self.masses[iso-1]
        denom = 2.0 * mass * self.rgrid2
        lam = self.state1._lambda
        sigma = self.state1.sigma
        omega = self.state1.omega
        spin = self.state1.spin
        num = (jrotn * (jrotn + 1.0) + spin * (spin + 1.0)) - (omega**2) - \
              (sigma**2) - (lam**2) + self.rot_correction

        # diag_values = self.Gy2 * (self.ugrid + (num / denom))
        diag_values = self.ygrid + (num / denom)
        np.fill_diagonal(self.matrix, diag_values)

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):
        """Calculate the radial function on a grid of points

        Args:
            ypar (array): The list of radial parameters
        """

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_pec_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_pec_on_grid(ypar)

        if self.model == 'morse':
            self._calculate_morse_pec_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_pec_on_grid(ypar)

        if self.model == 'emo':
            self._calculate_emo_pec_on_grid(ypar)

    def _calculate_pointwise_pec_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # TODO: The parameters in ypar should be passed in cm-1 units!!!
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_pec_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk_func = cs(self.rgrid, return_deriv=True)
        # stp, enp = (ch-1)*npnts, ch*npnts
        self.sk_grid[self.sind:self.eind, :] = sk_func

    def _calculate_morse_pec_on_grid(self, ypar):

        ypnts = ypar[self.sind:self.eind]
        # convert to au
        ypnts[0] = ypnts[0] * self.yunits[0]
        ypnts[1] = ypnts[1] * self.yunits[1]
        ypnts[2] = ypnts[2] * _utils.C_bohr
        ypnts[3] = ypnts[3] / _utils.C_bohr
        self.ygrid = self._morse_potential(ypnts)

    def _morse_potential(self, params):

        term_enr, diss_enr, aexp, req = params
        pwr = np.power((1.0 - np.exp(-aexp*(self.rgrid-req))), 2.0)
        morse_func = term_enr + (diss_enr * pwr)
        return morse_func

    def _calculate_custom_pec_on_grid(self, ypar):

        # TODO: check units
        ypnts = ypar[self.sind:self.eind] * self.yunits
        self.ygrid = self.custom_func(ypnts, self.rgrid)

    def _calculate_emo_pec_on_grid(self, ypar):

        ypnts = ypar[self.sind:self.eind]
        # TODO: convert to au
        self.ygrid = self._emo_potential(ypnts)

    def _emo_potential(self, params):

        term_enr, diss_enr, pexp, req = params[:4]
        bparams = np.array(params[4:])[::-1]

        yr_expr = self._calculate_emo_power_expression(req, pexp)

        bemo = np.polyval(bparams, yr_expr)
        pwr = np.power((1.0 - np.exp((-1.0*bemo)*(self.rgrid-req))), 2.0)
        vemo = term_enr + (diss_enr * pwr)

        return vemo

    def _calculate_mlr_pec_on_grid(self, ypar):

        ypnts = ypar[self.sind:self.eind]
        # TODO: convert to au
        self.ygrid = self._morse_long_range_potnetial(ypnts)

    def _morse_long_range_potnetial(self, params):

        nic, nbc, ncc, ndc = params[0:4]
        term_enr, diss_enr, pexp, qexp, rref, req, binf = params[:nic]

        bparams = np.array(params[nic:nic+nbc])[::-1]
        cparams = np.array(params[nic+nbc:nic+nbc+ncc])[::-1]
        dparams = np.array(params[nic+nbc+ncc:nic+nbc+ncc+ndc])[::-1]

        yrp = self._calculate_emo_power_expression(rref, pexp)
        yrq = self._calculate_emo_power_expression(rref, qexp)

        bmlj = yrp * binf + (1.0 - yrp) * np.polyval(bparams, yrq)
        ulrr = self._long_range_function(self.rgrid, cparams, dparams)
        ulrre = self._long_range_function(req, cparams, dparams)
        ulr = ulrr / ulrre
        vmlj = term_enr + diss_enr * np.power(1.0 - ulr * np.exp((-1.0*bmlj)*yrp), 2.0)

        return vmlj

    def _calculate_emo_power_expression(self, rref, power):

        numer = np.power(self.rgrid, power) - rref**power
        denom = np.power(self.rgrid, power) + rref**power

        return numer / denom

    def _long_range_function(self, rdist, cparams, dparams):

        # TODO: rewrite using numpy
        ulr = 0
        for i in range(0, cparams.shape[0]):
            ulr += dparams[i] * (cparams[i] / np.power(rdist, i+1))
        return ulr
