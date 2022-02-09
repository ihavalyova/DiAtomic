import numpy as np
import scipy as sp
import identify

from math import sqrt as _sqrt
from scipy.interpolate import CubicSpline as _CubicSpline
from os.path import splitext as _splitext
from interpolator import CSpline
from utils import C_hartree
from utils import Utils as _utils

__all__ = ['Hamiltonian']


class Hamiltonian:
    """Construct the Hamiltonian operator, solve the time-independant coupled
    system of radial Shrodinger equations and build a list of energy levels
    """

    def __init__(self, grid, diatomic_data, channels, couplings=[],
                 eig_decomp='lapack', lapack_driver='evr',
                 arpack_options=('LM', 6, None), is_weighted=False):
        """Initilize the Hamiltonian object

        Args:
            grid (object): The Grid object containing the grid parameters
            diatomic_data (object): a DiatomicData object containing the
                parameters about the specfied molecule.
            channels (list): The list of the defined channel objects.
            couplings (list, optional) The list of the defined coupling
                objects. Defaults to [].
            eig_decomp (str, optional): Which package to be used for eigen
                decomposition, LAPACK or ARPACK. Defaults to 'lapack'.
            lapack_driver (str, optional) Defines which LAPACK driver should be
                used when eig_decomp is set to "lapack". Valid options are:
                "ev", "evd", "evr", "evx". syevr() procedure is called by
                default. Defaults to 'evr'.
            arpack_options (tuple, optional)
                A list of optional parameters about the arpack diagonalization.
                Defaults to ('LM', 6, None).
            is_weighted : (bool, optional). whether to apply Watson's weighting
                procedure. Defaults to False.
        """

        self.ngrid = grid.ngrid
        self.rmin = grid.rmin
        self.rmax = grid.rmax
        self.rgrid = grid.rgrid
        self.rgrid2 = np.square(self.rgrid)
        self.rstep = grid.rstep
        self.solver = grid.solver

        self.molecule = diatomic_data.molecule
        self.mol_name = ''.join(
            filter(lambda x: not x.isdigit(), self.molecule[0]))
        self.jqnumbers = diatomic_data.jqnumbers
        self.pars = sorted(diatomic_data.symmetry)
        self.masses = diatomic_data.reduced_masses or diatomic_data.masses
        self.niso = diatomic_data.niso
        self.refj = diatomic_data.refj
        self.refE = diatomic_data.refE
        self.exp_data = diatomic_data.exp_data
        self.exp_file = diatomic_data.exp_file
        self.wavens_data = diatomic_data.wavens_data

        self.channels = channels
        self.couplings = couplings
        self.nch = len(self.channels)
        self.ncp = len(self.couplings)

        # mapping arrays
        self.Fy = grid.Fy
        self.Gy = grid.Gy

        self.msize = self.nch * self.ngrid
        self.dd = np.diag_indices(self.ngrid)
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
        self.nmax_levels = 10000

        # energy ranges when extracting term values
        # self.uenergy_range = None
        # self.lenergy_range = None

        self.kin_enr = KineticEnergy(self)
        self.pot_enr = PotentialEnergy(self)
        self.interact = Interaction(self)

        # matrix with the spline S functions
        self.sk_grid = np.zeros((self.msize, nmax_params))

        # used in the fit for initilization of the parameters
        channel_pars = diatomic_data.get_channel_parameters(channels)
        self.tot_ppts = channel_pars[0].shape[0]

        coupling_pars = diatomic_data.get_coupling_parameters(couplings)
        self.tot_cpts = coupling_pars[0].shape[0]

        self.ypar_init = np.concatenate((channel_pars[0], coupling_pars[0]))
        self.yfixed_init = np.concatenate((channel_pars[1], coupling_pars[1]))
        self.xunits = np.concatenate((channel_pars[2], coupling_pars[2]))
        self.yunits = np.concatenate((channel_pars[3], coupling_pars[3]))

        # get other parameters and functions used in the fit
        self.unq_channel_inds = diatomic_data.unq_chind
        self.edit_channel_parameters = diatomic_data.edit_channel_parameters
        self.edit_coupling_parameters = diatomic_data.edit_coupling_parameters

        self.interpolate_functions(self.ypar_init)

        self.file_energies = f'energies_{self.mol_name}.dat'

    def interpolate_functions(self, ypar):
        """Interpolate coupling and channel functions on the grid

        Args:
            ypar (array): An array containing the channel/coupling parameters
        """

        self.pot_enr.calculate_channels_on_grid(ypar=ypar)
        self.ugrid = self.pot_enr.ugrid

        if self.ncp != 0:
            self.interact.calculate_couplings_on_grid(ypar=ypar)

        self.fgrid = self.interact.fgrid

    def get_energy_levels(self, ident_option=1, store_output=True,
                          out_file=None):

        self.ident_option = ident_option if ident_option in [0, 1] else 1

        if self.exp_data is not None:
            self.out_data = identify.identify_levels(self.calc_data,
                                                     self.exp_data,
                                                     self.nch,
                                                     ident_option)
            # self.save_full_energy_list(self.calc_data)
            # self.save_energy_list(stats)
        else:
            raise AttributeError(
                'Experimental energy terms should be provided.')

        return self.out_data

    def solve(self, energy_subset_index=None, energy_subset_value=None):
        """A general function which solves the coupled system of Shcrodinger equations

        Args:
            energy_subset_index (tuple, optional): Specifies a subset of
                energies by providing two indices: the initial and the final
                index of the coumputed eigenvalues. Defaults to None.
            energy_subset_value (tuple, optional): Specifies a subset of
                energies by providing two energy values: the start and the end
                value of the coumputed eigenvalues. Defaults to None.

        Returns:
            tuple: A tuple of the computed eigenvalues and eigenvectors
        """

        self.energy_subset_value = energy_subset_value
        self.energy_subset_index = energy_subset_index

        # count the total number of eigenvalues
        self.ecount = 0

        evalues_all = np.zeros((self.nmax_levels, 9+self.nch))
        evecs_all = np.zeros((self.msize, self.nmax_levels))

        for niso, iso in enumerate(self.niso):
            tmatrix = self.kin_enr.calculate_kinetic_energy(iso)

            for par in self.pars:
                shiftE = [0.0]  # pass by reference; make it self?
                for jrotn in self.jqnumbers:

                    # build and diagonalize the hamiltonian matrix
                    hmatrix = self.build_hamiltonian(iso, par, jrotn, tmatrix)
                    evalues, evecs = self.diagonilize[self.eig_decomp](hmatrix)

                    # arange and store the quantum numbers & labels for levels
                    nevalues = evalues.shape[0]
                    levels = self._arange_levels(
                        jrotn, par, iso, evalues, evecs, shiftE=shiftE)
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

        self.interpolate_functions(ypar)
        self.solve(energy_subset_index=self.energy_subset_index,
                   energy_subset_value=self.energy_subset_value)
        return self.get_energy_levels(self.ident_option)

    def extract_terms_in_range(self, uenergy=None, lenergy=None, usymm=None,
                               lsymm=None, uj=None, lj=None, ustate=None,
                               lstate=None):

        # TODO: it may happen that the number of levels within the specified
        # ranges will vary during the fit

        eind, jind, pind, sind = 1, 2, 3, -4

        # filter by energy
        self.uenergy = uenergy
        self.lenergy = lenergy
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
            tmatrix (array): the kinetic energy matrix

        Returns:
            array: the final constructed Hamiltonian matrix
        """
        vmatrix = self.pot_enr.calculate_potential_energy(jrotn, par, iso)
        imatrix = self.interact.calculate_coupling_matrix(jrotn, par, iso)
        hmatrix = tmatrix + vmatrix + imatrix
        # evalues, evecs = self.diagonilize[self.eig_decomp](hmatrix)

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
        if self.energy_subset_index:
            emini, emaxi = self.energy_subset_index[0]
            emaxi = self.energy_subset_index[-1]
            subset_index = (emini, emaxi-1)

        elif self.energy_subset_value:
            eminv = self.energy_subset_value[0] / C_hartree
            emaxv = self.energy_subset_value[-1] / C_hartree
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
            return_eigenvectors=True
        )

        return evalues, evecs

    def arange_levels(self, jrotn, par, iso, evalues, evecs, shiftE=[0.0]):
        # nevalues = evalues.shape[0]
        # levels = self._arange_levels(
        #   jrotn, par, iso, evalues, evecs, shiftE=shiftE)
        # self.calc_data[self.ecount:nevalues+self.ecount, :] = levels
        # self.evecs_matrix[:, self.ecount:nevalues+self.ecount] = evecs
        # self.ecount += nevalues
        pass

    def _arange_levels(self, jrotn, par, iso, evalues, evecs, shiftE=[0.0]):
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
            shiftE (list, optional): The shift energy. Defaults to [0.0].

        Returns:
            array: The energy levels with the corresponding information
            aranged by columns
        """

        ids = np.arange(self.ecount+1, evalues.shape[0]+self.ecount+1)

        if self.refj and jrotn == self.jqnumbers[0]:
            shiftE[0] = evalues[0]
        elif self.refE:
            shiftE[0] = self.refE

        evalues_shifted = (evalues - shiftE[0]) * C_hartree

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
            map(lambda x: self.channels[x-1].nlambda, np.int64(states)),
            dtype=np.int64)

    def _get_omega_values(self, states):

        return np.fromiter(
            map(lambda x: self.channels[x-1].omega, np.int64(states)),
            dtype=np.float64)

    def save_energy_list(self, stats=None):

        if stats is None:
            stats = _utils.calculate_stats(self.out_data[:, 8],
                                           self.out_data[:, 7],
                                           self.out_data[:, 10],
                                           self.is_weighted)
        # add id column
        outdata_num = np.c_[np.arange(1, self.out_data.shape[0]+1),
                            self.out_data]

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

        np.savetxt(self.energy_out_file, outdata_num, header=header,
                   footer=footer, fmt=fmt)

    def save_full_energy_list(self, calc_energies, filename=None):
        """Stores the complete list of computed energy levels in external file

        Args:
            calc_energies (array): the array containing the energy levels
            filename (str, optional): The name of the file. Defaults to None.
        """

        nrows = self.calc_data.shape[1]
        cols = [0, -3] + list(range(1, nrows-3)) + [-2, -1]
        # calc_data_out = self.calc_data[:, cols]
        calc_energies_out = calc_energies[:, cols]

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

        fname_def, fext_def = _splitext(self.file_energies)
        output_file = fname_def + '_predicted' + fext_def
        output_file = filename or output_file

        np.savetxt(output_file, calc_energies_out, header=header, fmt=fmt)

    def _sort_predicted_energies(self, cols=[]):

        out = [self.calc_data[:, col] for col in reversed(cols)]
        self.calc_data = self.calc_data[np.lexsort(tuple(out))]

    def sort_predicted_energies(self, cols=[]):

        self._sort_predicted_energies(cols)
        self.save_full_energy_list()

    def _sort_output_energies(self, cols=[]):

        cols = [i-2 for i in cols]
        out = [self.out_data[:, col] for col in reversed(cols)]
        self.out_data = self.out_data[np.lexsort(tuple(out))]

    def sort_output_energies(self, cols=[]):

        self._sort_output(cols)
        self.save_energy_list()

    def get_predicted_data(self):
        """Get the complete list of computed energy levels

        Returns:
            array: the computed energy levels
        """

        return self.calc_data

    def get_output_data(self):

        return self.out_data

    def get_hamiltonian(self):
        """Get the Hamiltonian matrix for the last computed J, symmetry
        and isotopologue

        Returns:
            array: The Hamiltonian matrix
        """

        return self.hamiltonian


class KineticEnergy:
    """Calculate the kinetic energy operator
    """
    def __init__(self, H):

        self.rmin = H.rmin
        self.rmax = H.rmax
        self.ngrid = H.ngrid
        self.nch = H.nch
        self.solver_method = H.solver
        self.masses = H.masses
        self.Fy = H.Fy

        self.T = np.zeros((H.msize, H.msize))

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

        if self.solver_method.lower() == 'sinc':
            return self._calculate_kinetic_energy_sinc(mass)

        elif self.solver_method.lower() == 'fourier':
            return self._calculate_kinetic_energy_fourier(mass)

        elif self.solver_method.lower() == 'fd5':
            return self._calculate_kinetic_energy_FD5(mass)

        else:
            raise ValueError(
                f'{self.solver_method} is not allowed solver method...')

    def _calculate_kinetic_energy_fourier(self, mass):
        """Calculate the kinetic energy operator using the Fourir basis

        Args:
            mass (double): The mass value

        Returns:
            array: The computed kinetic energy matrix
        """

        length = self.rmax - self.rmin
        n2 = (self.ngrid**2 + 2.0) / 6.0
        pc = 1.0
        h = 2.0 * np.pi * pc
        ml = 4.0 * mass * length**2
        h2 = h**2

        ij_grid = np.mgrid[:self.ngrid, :self.ngrid]
        ij_mtx = ij_grid[0] - ij_grid[1]  # TODO: check this
        di = np.diag_indices(self.ngrid)
        sinf_mtx = np.sin((ij_mtx * np.pi) / self.ngrid)

        # set diagonal to some temporary nonzero value
        sinf_mtx[di] = 1.0

        self.T[:self.ngrid, :self.ngrid] = \
            (np.power(-1.0, ij_mtx) * h2) / (ml * np.power(sinf_mtx, 2))

        self.T[di] = (h2 / ml) * n2

        for ch in range(2, self.nch+1):
            ind1 = (ch-1)*self.ngrid
            ind2 = ch*self.ngrid
            self.T[ind1:ind2, ind1:ind2] = self.T[:self.ngrid, :self.ngrid]

        return self.T

    def _calculate_kinetic_energy_sinc(self, mass):
        """Calculate the kinetic energy operator using the sinc basis

        Args:
            mass (double): The mass value

        Returns:
            array: The computed kinetic energy matrix
        """

        pc = 1.0
        pi2 = np.pi ** 2
        h = (self.rmax - self.rmin) / (self.ngrid-1)
        h2 = (pc ** 2) / (2.0 * mass)

        ij_grid = np.mgrid[:self.ngrid, :self.ngrid]
        ij_diff = ij_grid[0] - ij_grid[1]
        ij_sum = ij_grid[0] + ij_grid[1]
        di = np.diag_indices(self.ngrid)

        # set diagonal to some temporary nonzero value
        ij_diff[di] = 1.0

        self.T[:self.ngrid, :self.ngrid] = \
            (2.0 * np.power(-1.0, ij_sum)) / np.power(ij_diff, 2)
        self.T[di] = pi2 / 3.0
        self.T = self.T * (h2 / h**2)

        for ch in range(2, self.nch+1):
            ind3 = (ch-1)*self.ngrid
            ind4 = ch*self.ngrid
            self.T[ind3:ind4, ind3:ind4] = self.T[:self.ngrid, :self.ngrid]

        self.T[di] = self.T[di] - (h2 * self.Fy)

        return self.T

    def _calculate_kinetic_energy_FD5(self, mass):
        """Calculate the kinetic energy operator using the finite difference method

        Args:
            mass (double): The mass value

        Returns:
            array: The computed kinetic energy matrix
        """

        # the first and last 2 eigenvalues are wrong
        pc = 1.0
        h2 = (pc ** 2) / (2.0 * mass)
        step = (self.rmax - self.rmin) / (self.ngrid-1)

        d0 = np.empty(self.ngrid)
        d0.fill(5.0/2.0)

        d1 = np.empty(self.ngrid-1)
        d1.fill(-4.0/3.0)

        d2 = np.empty(self.ngrid-2)
        d2.fill(1.0/12.0)

        self.T[:self.ngrid, :self.ngrid] = \
            (h2/(step**2)) * self._five_diagonal_matrix(d0, d1, d2)

        corner_coef = 29.0 / 12.0
        self.T[0, 0] = corner_coef
        self.T[self.ngrid-1, self.ngrid-1] = corner_coef

        for ch in range(2, self.nch+1):
            ind3 = (ch-1)*self.ngrid
            ind4 = ch*self.ngrid
            self.T[ind3:ind4, ind3:ind4] = self.T[:self.ngrid, :self.ngrid]

        return self.T

    def _five_diagonal_matrix(self, a, b, c, k1=0, k2=1, k3=2):

        return np.diag(a, k1) + np.diag(b, k2) + \
               np.diag(b, -k2) + np.diag(c, k3) + \
               np.diag(c, -k3)


class PotentialEnergy:

    def __init__(self, H):

        self.nch = H.nch
        self.ncp = H.ncp
        self.ngrid = H.ngrid
        self.rgrid = H.rgrid
        self.rgrid2 = H.rgrid2
        self.channels = H.channels
        self.couplings = H.couplings
        self.masses = H.masses
        self.dd = H.dd
        self.Gy2 = np.power(H.Gy, 2)

        self.ugrid = np.zeros(H.msize)
        self.V = np.zeros((H.msize, H.msize))

    def calculate_potential_energy(self, jrotn, par, iso):

        mass = self.masses[iso-1]
        denom = 2.0 * mass * self.rgrid2

        for ch in range(0, self.nch):
            # ctype = self.channels[ch].model
            lam = self.channels[ch].nlambda
            sigma = self.channels[ch].nsigma
            omega = self.channels[ch].omega
            spin = (self.channels[ch].mult - 1.0) / 2.0
            rot_correction = self.channels[ch].rot_correction

            num = (jrotn * (jrotn + 1.0) + spin * (spin+1.0)) - \
                (omega**2) - (sigma**2) - (lam**2) + rot_correction

            i1 = ch*self.ngrid
            i2 = (ch+1)*self.ngrid

            self.V[i1:i2, i1:i2][self.dd] = \
                self.Gy2 * (self.ugrid[i1:i2] + (num / denom))

        return self.V

    def calculate_channels_on_grid(self, ypar=None):

        for ch in range(self.nch):
            # built-in python cubic spline function
            if self.channels[ch].model.lower() == 'pointwise':
                self._calculate_pointwise_pec_on_grid(ch, ypar)

            # custom implementation of cubic spline function
            if self.channels[ch].model.lower() == 'cspline':
                self._calculate_cspline_pec_on_grid(ch, ypar)

            # Morse potential function
            if self.channels[ch].model.lower() == 'morse':
                self._calculate_Morse_pec_on_grid(ch, ypar)

            # EMO potential function
            if self.channels[ch].model.upper() == 'emo':
                self._calculate_EMO_pec_on_grid(ch, ypar)

            # MLR potential function
            if self.channels[ch].model.upper() == 'mlr':
                self._calculate_MLR_pec_on_grid(ch, ypar)

            # custom (user-specified) potential function
            if self.channels[ch].model.lower() == 'custom':
                self._calculate_custom_pec_on_grid(ch, ypar)

    def _calculate_pointwise_pec_on_grid(self, ch, ypar):

        sind = self.channels[ch].start_index
        eind = self.channels[ch].end_index
        xpnts = self.channels[ch].rpoints
        ypnts = ypar[sind:eind]
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ugrid[ch*self.ngrid:(ch+1)*self.ngrid] = cs(self.rgrid)

    def _calculate_cspline_pec_on_grid(self, ch, ypar):

        sind = self.channels[ch].start_index
        eind = self.channels[ch].end_index
        xpnts = self.channels[ch-1].R
        npnts = self.channels[ch-1].npnts
        ypnts = ypar[sind:eind]

        cs = CSpline(xpnts, ypnts)
        index1, index2 = (ch-1)*self.ngrid, ch*self.ngrid
        self.ugrid[index1:index2], sk = cs(self.rgrid, return_deriv=True)

        stp, enp = (ch-1)*npnts, ch*npnts
        self.sk_grid[(ch-1)*self.ngrid:ch*self.ngrid, stp:enp] = sk

    def _calculate_Morse_pec_on_grid(self, ch, ypar, ci):

        ypnts = self.channels[ch-1].U
        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = self._VMorse(ypnts)

    def _VMorse(self, params):

        Te, De, a, re = params
        return Te + De*np.power((1.0 - np.exp(-a*(self.rgrid-re))), 2.0)

    def _calculate_EMO_pec_on_grid(self, ch, ypar, ci):

        ypnts = self.channels[ch-1].U
        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = self._VEMO(ypnts)

    def _VEMO(self, params):

        Te, De, p, re = params[:4]
        bparams = np.array(params[4:])[::-1]

        yr = self.calculate_emo_power_expression(re, p)

        bemo = np.polyval(bparams, yr)
        pwr = np.power((1.0 - np.exp((-1.0*bemo)*(self.rgrid-re))), 2.0)
        vemo = Te + De * pwr

        return vemo

    def _calculate_MLR_pec_on_grid(self, ch, ypar, ci):

        ni = self.channels[ch-1].ni
        nb = self.channels[ch-1].nb
        nc = self.channels[ch-1].nc
        nd = self.channels[ch-1].nd

        nall = (ni, nb, nc, nd)

        ypnts = self.channels[ch-1].U
        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = self._VMLR(ypnts, nall)

    def _VMLR(self, params, nparams):

        ni, nb, nc, nd = nparams
        Te, De, p, q, rref, re, binf = params[:ni]

        bparams = np.array(params[ni:ni+nb])[::-1]
        cparams = np.array(params[ni+nb:ni+nb+nc])[::-1]
        dparams = np.array(params[ni+nb+nc:ni+nb+nc+nd])[::-1]

        yrp = self.calculate_emo_power_expression(rref, p)
        yrq = self.calculate_emo_power_expression(rref, q)

        bmlj = yrp * binf + (1.0 - yrp) * np.polyval(bparams, yrq)
        ulrr = self._long_range_function(self.rgrid, cparams, dparams)
        ulrre = self._long_range_function(re, cparams, dparams)
        ulr = ulrr / ulrre
        vmlj = Te + De * np.power(1.0 - ulr * np.exp((-1.0*bmlj)*yrp), 2.0)

        return vmlj

    def calculate_emo_power_expression(self, rr, power):

        numer = np.power(self.rgrid, power) - rr**power
        denom = np.power(self.rgrid, power) + rr**power

        return numer / denom

    def _long_range_function(self, r, cparams, dparams):

        # TODO: rewrite using numpy
        ulr = 0
        for i in range(0, cparams.shape[0]):
            ulr += dparams[i] * (cparams[i] / np.power(r, i+1))
        return ulr

    def _calculate_custom_pec_on_grid(self, ch, ypar, pranges, index):

        ypnts = self.channels[ch-1].U

        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
            self.channels[ch-1].cfunc(ypnts, self.rgrid)


class Interaction:

    def __init__(self, H):

        self.ngrid = H.ngrid
        self.nch = H.nch
        self.ncp = H.ncp
        self.channels = H.channels
        self.couplings = H.couplings
        self.rgrid = H.rgrid
        self.rgrid2 = H.rgrid2
        self.masses = H.masses
        self.msize = H.msize
        self.Gy = H.Gy
        self.dd = H.dd
        self.msize = H.msize
        self.fgrid = np.zeros(self.ncp * self.ngrid)

        self.countl = 0
        self.cp_map = {}
        max_rotnj = 200

        # TODO: the upper limit is not correct
        for cp in range(0, self.ncp+6):
            self.cp_map[cp] = np.zeros((self.msize, max_rotnj))

    def calculate_coupling_matrix(self, jrotn, par, iso):

        self.jrotn = jrotn
        mass = self.masses[iso-1]
        jjrotn = self.jrotn*(self.jrotn + 1.0)

        pert_matrix = np.zeros((self.msize)**2).reshape(self.msize, self.msize)

        interact_keys = self.define_interaction_keys()

        for cp in range(0, len(self.couplings)):
            cprops = zip(self.couplings[cp].interact,
                         self.couplings[cp].coupling,
                         self.couplings[cp].multiplier)

            ycs = self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid]

            for countc, (inter, ctype, m) in enumerate(cprops):
                ch1, ch2 = inter[:2]
                args = self.get_quantum_numbers(ch1, ch2)
                cfunc = interact_keys[ctype](jjrotn, mass, m, par, args)  # Gy2
                row1, row2 = (ch1-1)*self.ngrid, ch1*self.ngrid
                col1, col2 = (ch2-1)*self.ngrid, ch2*self.ngrid

                pert_matrix[row1:row2, col1:col2][self.dd] += cfunc * ycs
                # self.cp_map[cp+countc][row1:row2, self.countl] = cfunc

        self.countl += 1

        return pert_matrix

    def get_couplings_on_grid(self):
        return self.fgrid.reshape(self.ncp, self.ngrid).T  # * C_hartree

    def calculate_couplings_on_grid(self, ypar=None):

        for cp in range(self.ncp):
            # built-in cubic spline function
            if self.couplings[cp].model.lower() == 'pointwise':
                self._calculate_pointwise_coupling_on_grid(cp, ypar)

            # custom implementation of cubic spline function
            if self.couplings[cp].model.lower() == 'cspline':
                self._calculate_cspline_coupling_on_grid(cp, ypar)

            # constant coupling value
            if self.couplings[cp].model.lower() == 'constant':
                self.calculate_constant_coupling_on_grid(cp, ypar)

            # custom coupling function
            if self.couplings[cp].model.lower() == 'custom':
                self._calculate_custom_coupling_on_grid(cp, ypar)

            # self._countp += self.couplings[cp].npnts

    def _calculate_pointwise_coupling_on_grid(self, cp, ypar):

        xpnts = self.couplings[cp].xc
        ypnts = \
            ypar[self.couplings[cp].start_index:self.couplings[cp].end_index]

        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, cp, ypar):

        xpnts = self.couplings[cp].xc
        ypnts = self.couplings[cp].yc

        cs = CSpline(xpnts, ypnts)
        index = cp*self.ngrid, (cp+1)*self.ngrid

        self.fgrid[index[0]:index[1]], sk = cs(self.rgrid, return_deriv=True)

    # def calculate_coupling_points(self, cp, yrange, ypar):
    #     yrange['end'] = yrange['start'] + self.couplings[cp].npnts
    #     ypnts = ypar[yrange['start']:yrange['end']]
    #     yrange['start'] += self.couplings[cp].npnts
    #     return ypnts

    def calculate_constant_coupling_on_grid(self, cp, ypar):
        # st = Channel.totp_unique
        # en = st + self.couplings[cp].npnts
        # ygrid = ypar[st:en]
        # st += self.couplings[cp].npnts
        # self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = ygrid
        pass

    def _calculate_custom_coupling_on_grid(self, cp, ypar):

        ypnts = self.couplings[cp].yc
        # self.calculate_coupling_points(cp, yrange, ypar)

        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = \
            self.couplings[cp].cfunc(ypnts, self.rgrid)

    def get_quantum_numbers(self, ch1, ch2):

        return {
            'lm1': self.channels[ch1-1].nlambda,
            'sg1': self.channels[ch1-1].nsigma,
            'om1': self.channels[ch1-1].omega,
            'lm2': self.channels[ch2-1].nlambda,
            'sg2': self.channels[ch2-1].nsigma,
            'om2': self.channels[ch2-1].omega,
            's1': (self.channels[ch1-1].mult-1)/2,
            's2': (self.channels[ch2-1].mult-1)/2
        }

    def define_interaction_keys(self):

        return {
            'spin-orbit': self.spin_orbit_interaction,
            'lj': self.lj_interaction,
            'lj_parity': self.lj_parity_interaction,
            'sj': self.sj_interaction,
            'sj_parity': self.sj_parity_interaction,
            'sl': self.spin_electornic_interaction,
            'lambdadf': self.lambda_doubling_f_parity,
            'lambdade': self.lambda_doubling_e_parity,
            'lambdad': self.lambda_doubling,
            'ndbobc': self.NDBOBC,
            'dbobc': self.DBOBC,
            'spin-rot': self.spin_rotation_interaction,
            'spin-spin': self.spin_spin_interaction,
        }

    def spin_orbit_interaction(self, jjrotn, mass, m, par, args):
        """Calculate diagonal and off-diagonal Spin-Orbit matrix element
            ``<State1| SO |State1>`` = multiplier * A(R)
            ``<State1| SO |State2>`` = multiplier * alpha(R)

           with selection rules:

        Args:
            jjrotn (float): parameter with value equal to J(J+1)
            mass (float): the mass value
            m (float): multiplier for the R-dependent functions
            par (int): specify e- or f-symmetry
            args (list): the electronic and rotational quantum numbers

        Returns:
            array: the computed matrix elements
        """

        socoef = m * (self.rule_SOdiag(args) or self.rule_SOnondiag(args))

        # return (ycs / C_hartree) * socoef
        return socoef

    def rule_SOdiag(self, args):

        if args['lm1'] == args['lm2'] and \
           args['sg1'] == args['sg2'] and \
           args['om1'] == args['om2']:
            return 1

        return 0

    def rule_SOnondiag(self, args):

        if args['lm1'] != args['lm2'] or \
           args['sg1'] != args['sg2'] or \
           args['om1'] != args['om2']:
            return 1

        return 0

    def lj_interaction(self, jjrotn, mass, m, par, args):
        """Calculate the matrix elements of LJ operator which are:
           ``<State1| LJ |State1>`` =
           ``<State1| LJ |State2>`` =

           with selection rules:

        Args:
            jjrotn (float): parameter with value equal to J(J+1)
            mass (float): the mass value
            m (float): multiplier for the R-dependent functions
            par (int): specify e- or f-symmetry
            args (list): the electronic and rotational quantum numbers

        Returns:
            array: the computed matrix elements
        """

        qexpression = 1.0

        rule_omegaj = self.rule_lj_interaction(args) and \
            _sqrt(jjrotn - (args['om1'] * args['om2']))

        ljcoef = m * qexpression * rule_omegaj

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = -1.0

        return sign * brot * ljcoef

    def rule_lj_interaction(self, args):

        # to check the abs values
        omega_diff = abs(args['om1'] - args['om2'])
        lambda_diff = abs(args['lm1'] - args['lm2'])

        rule_lj = \
            ((omega_diff == lambda_diff) == 1.0) or \
            ((omega_diff == lambda_diff) == -1.0)

        # rule_lj = \
        #   ((args['om1'] - args['om2']) == \
        #   (args['lm1'] - args['lm2']) == 1.0) or \
        #    ((args['om1'] - args['om2']) == \
        #   (args['lm1'] - args['lm2']) == -1.0)
        # print('LJ', rule_lj)

        # # to check!!!
        sg1 = args['sg1']
        if args['om1'] == 0 and args['sg1'] < 0:
            sg1 = abs(args['sg1'])

        sg2 = args['sg2']
        if args['om2'] == 0 and args['sg2'] < 0:
            sg2 = abs(args['sg2'])

        # if args['sg1'] == args['sg2'] and \
        # args['s1'] == args['s2'] and rule_lj and \

        if sg1 == sg2 and args['s1'] == args['s2'] and rule_lj and \
           (args['om1'] <= self.jrotn and args['om2'] <= self.jrotn+1.0):
            return 1

        return 0

    def lj_parity_interaction(self, jjrotn, mass, m, par, args):

        qexpression = _sqrt(jjrotn + args['om1'] * args['om2'])

        ljcoef = m * qexpression * self.rule_lj_parity(args)

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0  # f-parity
        if par == 1:
            sign = -1.0  # e-parity

        return sign * brot * ljcoef

    def rule_lj_parity(self, args):

        """
        This interaction is allowed only between Sigma-Pi or Sigma-Sigma states
        with even multiplicity with Lambda <= 1 and abs(Omega) = 0.5
        """

        rule_even_mult = (2 * args['sg1'] + 1) % 2 == 0 and \
            (2 * args['sg2'] + 1) % 2 == 0

        if 0 <= args['lm1'] <= 1 and 0 <= args['lm2'] <= 1 and \
           (not args['lm1'] == args['lm2'] == 1) and rule_even_mult and \
           args['om1'] == args['om2'] == 0.5:
            return 1

        return 0

    def sj_interaction(self, jjrotn, mass, m, par, args):
        """Calculate the matrix elements of SJ operator which are:
            ``<State1| SJ |State1>`` =
            ``<State1| SJ |State2>`` =

           with selection rules:

        Args:
            jjrotn (float): parameter with value equal to J(J+1)
            mass (float): the mass value
            m (float): multiplier for the R-dependent functions
            par (int): specify e- or f-symmetry
            args (list): the electronic and rotational quantum numbers

        Returns:
            array: the computed matrix elements
        """

        qexpression = _sqrt(
            args['s1'] * (args['s1'] + 1) - args['sg1'] * args['sg2']
        )

        # TODO: check this!
        rule_omegaj = \
            self.rule_sj_interaction(args) and \
            _sqrt(jjrotn - (args['om1'] * args['om2']))

        sjcoef = m * qexpression * rule_omegaj

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = -1.0

        return sign * brot * sjcoef

    def rule_sj_interaction(self, args):

        omega_diff = args['om1'] - args['om2']
        sigma_diff = args['sg1'] - args['sg2']
        rule_sj = \
            (omega_diff == sigma_diff == 1.0) or \
            (omega_diff == sigma_diff == -1.0)

        if args['lm1'] == args['lm2'] and args['s1'] == args['s2'] and \
           (args['om1'] <= self.jrotn and args['om2'] <= self.jrotn+1.0) and \
           rule_sj:
            return 1

        return 0

    def sj_parity_interaction(self, jjrotn, mass, m, par, args):
        """
            only between Sigma states
        """

        qexpression = _sqrt(jjrotn + args['om1'] * args['om2'])

        sjcoef = m * qexpression * self.rule_sj_parity(args)

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0  # f-parity
        if par == 1:
            sign = -1.0  # e-parity

        return sign * brot * sjcoef

    def rule_sj_parity(self, args):

        if args['lm1'] == args['lm2'] == 0 and args['s1'] == args['s2']:
            return 1

        return 0

    def spin_electornic_interaction(self, jjrotn, mass, m, par, args):
        """Calculate the matrix elements of SL operator which are:
            ``<State1| SL |State1>`` =
            ``<State1| SL |State2>`` =

           with selection rules:

        Args:
            jjrotn (float): parameter with value equal to J(J+1)
            mass (float): the mass value
            m (float): multiplier for the R-dependent functions
            par (int): specify e- or f-symmetry
            args (list): the electronic and rotational quantum numbers

        Returns:
            array: the computed matrix elements
        """

        qexpression = _sqrt(
            (args['s1'] * (args['s1'] + 1)) - args['sg1'] * args['sg2'])

        slcoef = m * qexpression * self.rule_spin_electronic(args)

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0

        return sign * brot * slcoef

    def rule_spin_electronic(self, args):

        lambda_diff = args['lm1'] - args['lm2']
        sigma_diff = args['sg2'] - args['sg1']

        rule_ls = \
            (lambda_diff == sigma_diff == 1.0) or \
            (lambda_diff == sigma_diff == -1.0)

        if args['om1'] == args['om2'] and args['s1'] == args['s2'] and rule_ls:
            return 1

        return 0

    def lambda_doubling(self, jjrotn, mass, m, par, args):

        ldcoef = jjrotn / (2.0 * mass * self.rgrid2)**2

        return m * ldcoef

    def lambda_doubling_e_parity(self, jjrotn, mass, m, par, args):

        if par == 1:
            return self.lambda_doubling(jjrotn, mass, m, par, args)

        return np.zeros(self.rgrid2.shape[0])

    def lambda_doubling_f_parity(self, jjrotn, mass, m, par, args):

        if par == 0:
            return self.lambda_doubling(jjrotn, mass, m, par, args)

        return np.zeros(self.rgrid2.shape[0])

    def NDBOBC(self, jjrotn, mass, m, par, args):

        bocoef = jjrotn / (2.0 * mass * self.rgrid2)

        return m * bocoef * C_hartree

    def DBOBC(self, jjrotn, mass, m, par, args):
        """The diagonal BOB correction

        Args:
            jjrotn (float): parameter with value equal to J(J+1)
            mass (float): the mass value
            m (float): multiplier for the R-dependent functions
            par (int): specify e- or f-symmetry
            args (list): the electronic and rotational quantum numbers

        Returns:
            array: the computed correction
        """

        # return (1.0 / C_hartree) * m
        # return m * ((mass - self.masses[0]) / mass)
        return m * ((self.masses[0] - mass) / mass)

    def spin_rotation_interaction(self, jjrotn, mass, m, par, args):

        if self.rule_spin_rot_diag(args):

            qexpression = args['sg1']**2 - args['s1']*(args['s2'] + 1.0)

        elif self.rule_spin_rot_nondiag(args):
            ss1 = args['s1'] * (args['s1'] + 1)
            qexpression = \
                _sqrt(jjrotn - (args['om1'] * args['om2'])) * \
                _sqrt(ss1 - args['sg1'] * args['sg2'])

        srcoef = m * qexpression
        sign = 1.0

        return sign * srcoef

    def rule_spin_rot_diag(self, args):

        if args['om1'] == args['om2'] and \
           args['s1'] == args['s2'] and \
           args['sg1'] == args['sg2'] and \
           args['lm1'] == args['lm2']:
            return 1

        return 0

    def rule_spin_rot_nondiag(self, args):

        self.rule_sj_interaction(args)

    def spin_spin_interaction(self, jjrotn, mass, m, par, args):

        if self.rule_spin_spin_diag(args):
            qexpression = 3 * args['sg1']**2 - args['s1']*(args['s2'] + 1.0)

        elif self.rule_spin_spin_nondiag(args):
            qexpression = 1.0

        sscoef = m * qexpression
        sign = 1.0

        return sign * sscoef

    def rule_spin_spin_nondiag(self, args):

        lambda_diff = args['lm1'] - args['lm2']
        sigma_diff = args['sg2'] - args['sg1']

        # Delta Sigma = +/- 1
        rule_ss1_1so = \
            (lambda_diff == -1.0 * sigma_diff == 1.0) or \
            (lambda_diff == -1.0 * sigma_diff == -1.0)

        # Delta Sigma = +/- 2
        rule_ss2_2so = \
            (lambda_diff == -1.0 * sigma_diff == 2.0) or \
            (lambda_diff == -1.0 * sigma_diff == -2.0)

        rule_ss1 = \
            (args['s1'] - args['s2']) == 1.0 or \
            (args['s1'] - args['s2']) == -1.0 and \
            rule_ss1_1so

        rule_ss2 = \
            (args['s1'] - args['s2']) == 2.0 or \
            (args['s1'] - args['s2']) == -2.0 and \
            rule_ss2_2so

        if rule_ss1 or rule_ss2:
            return 1

        return 0

    def rule_spin_spin_diag(self, args):

        if args['om1'] == args['om2'] and \
           args['s1'] == args['s2'] and \
           args['sg1'] == args['sg2']:
            return 1

        return 0
