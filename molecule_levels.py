import math
import os

import numpy as np
import scipy as sp

from scipy.interpolate import CubicSpline
from interaction import Interaction
from molecule_data import Channel
from molecule_data import Coupling
from constants import Const
from utils import Utils
from grids import CSpline


class MoleculeLevels:

    def __init__(self, md, grid, channels, couplings=[], data=None):

        self.jqnumbers = md.jqnumbers
        self.pars = md.parities
        self.masses = md.reduced_masses or md.masses
        self.nisotopes = md.nisotopes
        self.refj = md.referencej
        self.refE = md.refE
        self.exp_file = md.exp_file
        self.ngrid = grid.ngrid
        self.rmin = grid.rmin
        self.rmax = grid.rmax
        self.rgrid = grid.rgrid
        self.solver = grid.solver
        self.Fy = grid.Fy
        # to multiply all R dependant functions with this array
        self.Gy2 = np.power(grid.Gy, 2)

        self.channels = channels
        self.nch = len(self.channels)

        self.couplings = couplings
        self.ncp = len(self.couplings)

        self.evalues_subseti = None
        self.evalues_subsetv = None

        self.eig_funcs = {
            'lapack': self.lapack_eig_decomposition,
            'arpack': self.arpack_eig_decomposition
        }

        self.exp_data = data
        self.ninit = 1

        self.rgrid2 = np.square(self.rgrid)
        self.ugrid = np.zeros(self.nch * self.ngrid)
        self.fgrid = np.zeros(self.ncp * self.ngrid)

        self.pot_enr = np.zeros((self.nch*self.ngrid, self.nch*self.ngrid))
        self.kin_enr = np.zeros_like(self.pot_enr)

        self.evals_predicted = np.array([], ndmin=2, dtype=np.float64)

        # total number of potential parameters
        self._ctot = 0

        # spline derivatives
        self.sderiv = None

        self.evals_file = 'eigenvalues.dat'
        self.predicted_file = 'evalues_predicted.dat'
        self.info_outfile = 'data_info.dat'

        self.evec_dir = os.path.join(Utils.get_current_dir(), 'eigenvectors')
        self.wavef_dir = os.path.join(Utils.get_current_dir(), 'wavefunctions')

        # self.evectors = np.zeros(
        # ((self.nch*self.ngrid)**2, len(jqnumbers)*self.vmax))

    def calculate_channels_on_grid(self, ypar=None):

        self._ctot = 0
        unique = set()
        pot_ranges = dict()

        # pass by reference
        index = [-1]

        for ch in range(1, self.nch+1):

            # built-in cubic spline function
            if self.channels[ch-1].model.lower() == Channel.models[1]:
                self.calculate_pointwise_pec_on_grid(
                    ch, ypar, pot_ranges, index
                )

            # cubic spline function using with custom implementation
            if self.channels[ch-1].model.lower() == Channel.models[2]:
                self.calculate_custom_pointwise_pec_on_grid(
                    ch, ypar, pot_ranges, index
                )

            # Morse potential function
            if self.channels[ch-1].model.lower() == Channel.models[3]:
                self.calculate_Morse_pec_on_grid(ch, ypar, pot_ranges, index)

            # EMO potential function
            if self.channels[ch-1].model.upper() == Channel.models[4]:
                self.calculate_EMO_pec_on_grid(ch, ypar, pot_ranges, index)

            # MLR potential function
            if self.channels[ch-1].model.upper() == Channel.models[5]:
                self.calculate_MLR_pec_on_grid(ch, ypar, pot_ranges, index)

            # Custom potential function
            if self.channels[ch-1].model.lower() == Channel.models[6]:
                self.calculate_custom_pec_on_grid(ch, ypar, pot_ranges, index)

            # determine _ctot value
            if self.channels[ch-1].filep not in unique:
                unique.add(self.channels[ch-1].filep)

                self._ctot += self.channels[ch-1].npnts

    def calculate_pointwise_pec_on_grid(self, ch, ypar, pot_ranges, index):

        # xpnts, ypnts = self.calculate_pec_points(ch, ypar, pot_ranges, index)

        xpnts = self.channels[ch-1].R

        # ypar contains parameters from all potentials and they
        # have to be transformed to the different channels
        if ypar is not None:
            pname = self.channels[ch-1].filep
            npt = self.channels[ch-1].npnts

            if pname in pot_ranges:
                st = pot_ranges[pname][0]
                en = pot_ranges[pname][1]
            else:
                st = index[0] + 1
                en = st + npt
                pot_ranges[pname] = (st, en)
                index[0] += npt

            ypnts = ypar[st:en]
        else:
            ypnts = self.channels[ch-1].U

        cs = CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = cs(self.rgrid)

    def calculate_custom_pointwise_pec_on_grid(self, ch, ypar, pot_ranges,
                                               index):

        xpnts = self.channels[ch-1].R
        ypnts = self.calculate_pec_points(ch, ypar, pot_ranges, index)

        cs = CSpline(xpnts, ypnts)
        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid], Sk = \
            cs(self.rgrid, return_deriv=True)

        self.sderiv = Sk * Const.hartree  # ???

    def calculate_pec_points(self, ch, ypar, pot_ranges, index):

        # map the potential parameters in ypar to each channel
        if ypar is not None:
            pname = self.channels[ch-1].filep
            npt = self.channels[ch-1].npnts

            if pname in pot_ranges:
                st = pot_ranges[pname][0]
                en = pot_ranges[pname][1]
            else:
                st = index[0] + 1
                en = st + npt
                pot_ranges[pname] = (st, en)
                index[0] += npt

            ypnts = ypar[st:en]
        else:
            ypnts = self.channels[ch-1].U

        return ypnts

    def calculate_Morse_pec_on_grid(self, ch, ypar, pot_ranges, index):

        ypnts = self.calculate_pec_points(ch, ypar, pot_ranges, index)

        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
            self.VMorse(ypnts, self.rgrid)

    def VMorse(self, params, rgrid):

        Te, De, a, re = params

        return Te + De*np.power((1.0 - np.exp(-a*(rgrid-re))), 2.0)

    def calculate_EMO_pec_on_grid(self, ch, ypar, pot_ranges, index):

        if ypar is not None:
            pass
        else:
            self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
                self.VEMO(self.channels[ch-1].U, self.rgrid)

    def VEMO(self, params, rgrid):

        Te, De, p, re = params[:4]
        bparams = np.array(params[4:])[::-1]

        yr = (np.power(rgrid, p) - re ** p) / (np.power(rgrid, p) + re ** p)

        bemo = np.polyval(bparams, yr)
        vemo = Te + De * np.power((1.0 - np.exp((-1.0*bemo)*(rgrid-re))), 2.0)

        return vemo

    def calculate_MLR_pec_on_grid(self, ch, ypar, pot_ranges, index):

        if ypar is not None:
            pass
        else:
            ni = self.channels[ch-1].ni
            nb = self.channels[ch-1].nb
            nc = self.channels[ch-1].nc
            nd = self.channels[ch-1].nd

            self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
                self.VMLR(self.channels[ch-1].U, self.rgrid, (ni, nb, nc, nd))

    def VMLR(self, params, rgrid, nparams):

        ni, nb, nc, nd = nparams

        Te, De, p, q, rref, re, binf = params[:ni]

        bparams = np.array(params[ni:ni+nb])[::-1]
        cparams = np.array(params[ni+nb:ni+nb+nc])[::-1]
        dparams = np.array(params[ni+nb+nc:ni+nb+nc+nd])[::-1]

        yrp = (np.power(rgrid, p) - rref**p) / (np.power(rgrid, p) + rref**p)
        yrq = (np.power(rgrid, q) - rref**q) / (np.power(rgrid, q) + rref**q)

        bmlj = yrp * binf + (1.0 - yrp) * np.polyval(bparams, yrq)

        ulrr = self._long_range_function(rgrid, cparams, dparams)
        ulrre = self._long_range_function(re, cparams, dparams)

        ulr = ulrr / ulrre

        vmlj = Te + De * np.power(1.0 - ulr * np.exp((-1.0*bmlj)*yrp), 2.0)

        return vmlj

    def _long_range_function(self, r, cparams, dparams):

        # TODO: Rewrite using numpy
        ulr = 0
        for i in range(0, cparams.shape[0]):
            ulr += dparams[i] * (cparams[i] / np.power(r, i+1))

        return ulr

    def calculate_custom_pec_on_grid(self, ch, ypar, pot_ranges, index):

        ypnts = self.calculate_pec_points(ch, ypar, pot_ranges, index)

        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
            self.channels[ch-1].cfunc(ypnts, self.rgrid)

    def calculate_couplings_on_grid(self, ypar=None):

        # pass by reference
        yrange = {'start': self._ctot, 'end': self._ctot}

        for cp in range(0, self.ncp):

            # built-in cubic spline function
            if self.couplings[cp].model.lower() == Coupling.models[1]:
                self.calculate_pointwise_couplings_on_grid(cp, yrange, ypar)

            # custom implementation of cubic spline function
            if self.couplings[cp].model.lower() == Coupling.models[2]:
                self.calculate_custom_pointwise_couplings_on_grid(
                    cp, yrange, ypar
                )

            # constant coupling value
            if self.couplings[cp].model.lower() == Coupling.models[3]:
                self.calculate_constant_coupling_on_grid(cp, ypar)

            # custom coupling function
            if self.couplings[cp].model.lower() == Coupling.models[4]:
                self.calculate_custom_coupling_on_grid(cp, yrange, ypar)

    def calculate_pointwise_couplings_on_grid(self, cp, yrange, ypar):

        xpnts = self.couplings[cp].xc
        ypnts = self.calculate_coupling_points(cp, yrange, ypar)

        cs = CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)

        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = cs(self.rgrid)

    def calculate_custom_pointwise_couplings_on_grid(self, cp, yrange, ypar):

        xpnts = self.couplings[cp].xc
        ypnts = self.calculate_coupling_points(cp, yrange, ypar)

        cs = CSpline(xpnts, ypnts)

        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = cs(self.rgrid)

    def calculate_coupling_points(self, cp, yrange, ypar):

        yrange['end'] = yrange['start'] + self.couplings[cp].npnts

        if ypar is None:
            ypnts = self.couplings[cp].yc
        else:
            ypnts = ypar[yrange['start']:yrange['end']]

        yrange['start'] += self.couplings[cp].npnts

        return ypnts

    def calculate_constant_coupling_on_grid(self, cp, ypar):

        ygrid = np.empty(self.rgrid.shape[0])

        st = self._ctot
        en = st + self.couplings[cp].npnts

        if ypar is None:
            ygrid.fill(self.couplings[cp].yc[0])
        else:
            ygrid = ypar[st:en]

        st += self.couplings[cp].npnts

        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = ygrid

    def calculate_custom_coupling_on_grid(self, cp, yrange, ypar):

        ypnts = self.calculate_coupling_points(cp, yrange, ypar)

        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = \
            self.couplings[cp].cfunc(ypnts, self.rgrid)

    def calculate_kinetic_energy(self, mass):

        if self.solver == 'sinc':
            self.calculate_kinetic_energy_sinc(mass)

        elif self.solver.lower() == 'fd5':
            self.calculate_kinetic_energy_FD5(mass)

        elif self.solver.lower() == 'fourier':
            self.calculate_kinetic_energy_fourier(mass)

        else:
            raise SystemExit(f'Error: Invalid solver {self.solver}')

    def lapack_eig_decomposition(self):

        evalues, evecs = sp.linalg.eigh(
            self.hmatrix,
            eigvals_only=False,
            overwrite_a=True,
            lower=False,
            subset_by_index=self.evalues_subseti,
            subset_by_value=self.evalues_subsetv,
            driver=self.lapack_driver.lower(),
            check_finite=False
        )

        return evalues, evecs

    def arpack_eig_decomposition(self):

        """ARPACK procedure for diagonalization of the Hamiltonian matrix

        Args:
            hmatrix (2D array): the Hamiltonian matrix

        Returns:
            tuple of numpy arrays: the computed eigenvalues and eigenvectors

        Remarks:
            ARAPCK procedure is the most efficient and suitable for finding
            the largest eigenvalues of a sparse matrix. If the smallest
            eigenvalues are desired then it is recommended to use a
            shift-invert mode. It transforms the eigenvalue problem to
            an eqivalent problem with shifted eigenvalues in which the
            small eigenvalues u become the large eigenvalues v: v = 1 / u
        """

        evalues, evecs = sp.sparse.linalg.eigsh(
            self.hmatrix,
            k=self.arpack_k,
            which=self.arpack_which.upper(),
            sigma=self.arpack_sigma,
            return_eigenvectors=True
        )

        return evalues, evecs

    def define_keywords(self):

        return {
            1: 'eig_decomp',
            2: 'energy_subset_index',
            3: 'energy_subset_value',
            4: 'identify',
            5: 'store_predicted',
            6: 'store_info',
            7: 'store_evecs',
            8: 'weighted',
            9: 'sort_output',
            10: 'sort_predicted',
            11: 'lapack_driver',
            12: 'arpack_which',
            13: 'arpack_k',
            14: 'arpack_sigma',
            15: 'file_evals_predicted',
            16: 'file_eigenvalues',
        }

    def calculate_levels(self, ypar=None, **kwargs):

        # ypar is only for internal use

        keys = self.define_keywords()

        self.eig_decomp = kwargs.get(keys[1]) or 'lapack'
        self.energy_subset_index = kwargs.get(keys[2])
        self.energy_subset_value = kwargs.get(keys[3])
        self.identify = kwargs.get(keys[4]) or 0
        self.store_predicted = kwargs.get(keys[5]) or False
        self.store_info = kwargs.get(keys[6]) or False
        self.store_evecs = kwargs.get(keys[7]) or False
        self.is_weighted = kwargs.get(keys[8]) or False
        self.sort_output = kwargs.get(keys[9])
        self.sort_predicted = kwargs.get(keys[10])
        self.file_predicted = kwargs.get(keys[15]) or self.predicted_file
        self.file_evals = kwargs.get(keys[16]) or self.evals_file

        # LAPACK options
        # the lapack procedure syevr() is used by default
        self.lapack_driver = kwargs.get(keys[11]) or 'evr'

        # ARPACK options
        self.arpack_which = kwargs.get(keys[12]) or 'LM'

        # in the arpack procedure k=6 by default
        self.arpack_k = kwargs.get(keys[13]) or 6

        # if sigma is set shift-invert mode is requested
        self.arpack_sigma = kwargs.get(keys[14])

        if self.energy_subset_index is not None:
            self.emini = self.energy_subset_index[0]
            self.emaxi = self.energy_subset_index[-1]
            self.evalues_subseti = (self.emini, self.emaxi-1)

        if self.energy_subset_value is not None:
            self.eminv = self.energy_subset_value[0] / Const.hartree
            self.emaxv = self.energy_subset_value[-1] / Const.hartree
            self.evalues_subsetv = (self.eminv, self.emaxv)

        chi2, rms, rmsd = self.calculate_eigenvalues(ypar)

        print(f'Chi Square = {chi2:<18.8f} | RMS = {rms:<15.8f}cm-1')
        # f'\nRMSD = {rmsd:.8f}'

    def calculate_eigenvalues(self, ypar=None):

        dd = np.diag_indices(self.ngrid)
        self.ninit = 1

        self.calculate_channels_on_grid(ypar=ypar)

        if self.ncp != 0:
            self.calculate_couplings_on_grid(ypar=ypar)

        interact = Interaction(self.ngrid, self.nch, self.rgrid2)
        count = 0

        for iso in self.nisotopes:
            self.evals_predicted = np.array([], ndmin=2, dtype=np.float64)

            mass = self.masses[iso-1]

            self.calculate_kinetic_energy(mass)

            for par in self.pars:
                # pass by reference
                shiftE = [0.0]

                for jrot in self.jqnumbers:

                    self.calculate_potential_energy(jrot, par, mass, dd)

                    self.hmatrix = self.kin_enr + self.pot_enr

                    pert_matrix = interact.get_perturbations(
                        jrot, mass, par, self.channels, self.couplings,
                        self.fgrid, dd, self.Gy2
                    )

                    # pmatrix = pmatrix*np.tile(self.Gy2,(pmatrix.shape[0],5))

                    self.hmatrix += pert_matrix

                    evalues, evecs = self.eig_funcs[self.eig_decomp]()

                    ccoefs = self.get_coupling_coefficients(evecs)
                    states = np.argmax(ccoefs, axis=1) + 1
                    vibnums = self.assign_vibrational_number(states)

                    if self.store_evecs:
                        # jrange = np.linspace(js, je,
                        # num=(je-js)+1.0, dtype=np.float, endpoint=True)

                        # self._save_eigenvectors(
                        #   jrot, par, iso,
                        #   evecs[:, 0:evalues.shape[0]],
                        #   vibnums, states
                        # )

                        os.makedirs(self.evec_dir, exist_ok=True)
                        fname = f'evector_{str(count)}.dat'
                        np.savetxt(
                            os.path.join(self.evec_dir, fname),
                            evecs, newline='\n', fmt='%15.8e'
                        )
                        count += 1

                    calc_data = self.arange_levels(
                        jrot, par, iso, shiftE, evalues
                    )

                    calc_data = np.column_stack(
                        (calc_data, ccoefs, states, vibnums)
                    ).reshape(evalues.shape[0], -1)

                    self.evals_predicted = np.append(
                        self.evals_predicted, calc_data
                    ).reshape(-1, calc_data.shape[1])

            nlambda = self.get_lambda_values()
            omega = self.get_omega_values()

            self.evals_predicted = np.column_stack(
                (self.evals_predicted, nlambda, omega)
            )

            chi2, rms, rmsd = 0.0, 0.0, 0.0

            if self.exp_data is not None:
                self.out_data = np.array([], dtype=np.float64)

                self.identify_levels()

                chi2, rms, rmsd = self.calculate_stats(
                    self.out_data[:, 8], self.out_data[:, 7],
                    self.out_data[:, 10], self.is_weighted
                )

                append = False if iso == 1 else True
                if self.sort_output:
                    self.sort_output_data()

                self.save_output_data(append, stats=(chi2, rms, rmsd))

                self.ninit += self.out_data.shape[0]

                if self.store_predicted:

                    if self.sort_predicted is not None:
                        self.sort_predicted_data()

                    self.save_all_calculated_data()

                if self.store_info:
                    self.save_additional_data()

                self.evecs_index = self.get_evecs_indicies()
            else:
                if self.sort_predicted is not None:
                    self.sort_predicted_data()

                self.save_all_calculated_data()

        return chi2, rms, rmsd

    def get_evecs_indicies(self):

        """find the indices of the computed eigenvectors corresponding
        to the selected eiegenvalues which will be used in the computation
        of partial derivatives based on the Hellman-Feynman theorem

        Returns:
            list: the indicies of the eigenvectors
        """
        # TODO: will not work for different isotopes

        evecs_index, _ = self.get_indices_of_matching_rows_simple(
            self.evals_predicted[:, [0, 1, 2, -3, -4, -1, -2]],
            self.out_data[:, [8, 2, 6, 1, -1, 3, 5]]
        )

        # all eigenvectors should be calculated
        evecs_index = evecs_index % self.ngrid
        evecs_index1 = self.group_increasing_sequences(evecs_index)

        evecs_index1 = list(evecs_index1)

        return evecs_index1

    def sort_predicted_data(self):

        # sort the data by 6 specified column numbers
        col1, col2, col3, col4, col5, col6 = self.sort_predicted

        self.evals_predicted = self.evals_predicted[np.lexsort((
                self.evals_predicted[:, col6], self.evals_predicted[:, col5],
                self.evals_predicted[:, col4], self.evals_predicted[:, col3],
                self.evals_predicted[:, col2], self.evals_predicted[:, col1]
        ))]

    def sort_output_data(self):

        # sort the data by 6 specified column numbers
        col1, col2, col3, col4, col5, col6 = self.sort_output

        self.out_data = self.out_data[np.lexsort((
                self.out_data[:, col6], self.out_data[:, col5],
                self.out_data[:, col4], self.out_data[:, col3],
                self.out_data[:, col2], self.out_data[:, col1]
        ))]

    def group_increasing_sequences(self, iseq):

        return [
            list(arr) for arr in np.split(
                iseq, np.where(np.diff(iseq) != 1)[0] + 1
            ) if arr.size > 0
        ]

    def get_lambda_values(self):

        return np.fromiter(map(
            lambda x: self.channels[x-1].nlambda,
            np.int64(self.evals_predicted[:, -2])),
            dtype=np.int64
        )

    def get_omega_values(self):

        return np.fromiter(map(
            lambda x: self.channels[x-1].omega,
            np.int64(self.evals_predicted[:, -2])),
            dtype=np.float64
        )

    def arange_levels(self, jrot, par, iso, shiftE, evalues):

        if jrot == self.jqnumbers[0] and self.refj is not None:
            shiftE[0] = evalues[0]

        elif self.refE is not None:
            shiftE[0] = self.refE

        evalues = (evalues - shiftE[0]) * Const.hartree

        calc_data = np.column_stack((
                evalues[:evalues.shape[0]],
                np.full((evalues.shape[0],), jrot),
                np.full((evalues.shape[0],), par),
                np.full((evalues.shape[0],), iso),
        ))

        return calc_data

    def calculate_kinetic_energy_fourier(self, mass):

        length = self.rmax - self.rmin
        n2 = (self.ngrid**2 + 2.0) / 6.0

        pc = 1.0
        h = 2.0 * np.pi * pc
        ml = 4.0 * mass * length**2
        h2 = h**2

        ij_grid = np.mgrid[0:self.ngrid, 0:self.ngrid]
        ij_mtx = ij_grid[0] - ij_grid[1]

        di = np.diag_indices(self.ngrid)

        sinf_mtx = np.sin((ij_mtx * np.pi) / self.ngrid)
        sinf_mtx[di] = 1.0  # set diagonal to some temporary nonzero value

        self.kin_enr[0:self.ngrid, 0:self.ngrid] = \
            (np.power(-1.0, ij_mtx) * h2) / (ml * np.power(sinf_mtx, 2))
        self.kin_enr[di] = (h2 / ml) * n2

        for ch in range(2, self.nch+1):
            ind1 = (ch-1)*self.ngrid
            ind2 = ch*self.ngrid

            self.kin_enr[ind1:ind2, ind1:ind2] = \
                self.kin_enr[0:self.ngrid, 0:self.ngrid]

    def calculate_kinetic_energy_sinc(self, mass):

        pc = 1.0
        pi2 = np.pi ** 2
        h = (self.rmax - self.rmin) / (self.ngrid-1)
        h2 = (pc ** 2) / (2.0 * mass)

        ij_grid = np.mgrid[0:self.ngrid, 0:self.ngrid]
        ij_diff = ij_grid[0] - ij_grid[1]
        ij_sum = ij_grid[0] + ij_grid[1]
        di = np.diag_indices(self.ngrid)

        # set diagonal to some temporary nonzero value
        ij_diff[di] = 1.0

        self.kin_enr[0:self.ngrid, 0:self.ngrid] = \
            (2.0 * np.power(-1.0, ij_sum)) / np.power(ij_diff, 2)
        self.kin_enr[di] = pi2 / 3.0
        self.kin_enr = self.kin_enr * (h2 / h**2)

        for ch in range(2, self.nch+1):
            ind3 = (ch-1)*self.ngrid
            ind4 = ch*self.ngrid

            self.kin_enr[ind3:ind4, ind3:ind4] = \
                self.kin_enr[0:self.ngrid, 0:self.ngrid]

        self.kin_enr[di] = self.kin_enr[di] - (h2 * self.Fy)

    def calculate_kinetic_energy_FD5(self, mass):

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

        self.kin_enr[0:self.ngrid, 0:self.ngrid] = \
            (h2/(step**2)) * self.form_five_diagonal_symm_matrix(d0, d1, d2)

        corner_coef = 29.0 / 12.0
        self.kin_enr[0, 0] = corner_coef
        self.kin_enr[self.ngrid-1, self.ngrid-1] = corner_coef

        for ch in range(2, self.nch+1):
            ind3 = (ch-1)*self.ngrid
            ind4 = ch*self.ngrid

            self.kin_enr[ind3:ind4, ind3:ind4] = \
                self.kin_enr[0:self.ngrid, 0:self.ngrid]

    def form_five_diagonal_symm_matrix(self, a, b, c, k1=0, k2=1, k3=2):

        return np.diag(a, k1) + np.diag(b, k2) + np.diag(b, -k2) + \
                np.diag(c, k3) + np.diag(c, -k3)

    def calculate_potential_energy(self, jrot, par, mass, dd):

        denom = 2.0 * mass * self.rgrid2

        for ch in range(0, self.nch):
            # ctype = self.channels[ch].model
            lam = self.channels[ch].nlambda
            sigma = self.channels[ch].nsigma
            omega = self.channels[ch].omega
            spin = (self.channels[ch].mult - 1.0) / 2.0
            rot_correction = self.channels[ch].rot_correction

            num = (jrot * (jrot + 1.0) + spin * (spin+1.0)) - \
                (omega**2) - (sigma**2) - (lam**2) + rot_correction

            i1 = ch*self.ngrid
            i2 = (ch+1)*self.ngrid

            self.pot_enr[i1:i2, i1:i2][dd] = \
                self.Gy2 * (self.ugrid[i1:i2] + (num / denom))

    def get_coupling_coefficients(self, evecs):
        # used to assign state to the computed levels
        # based on the largest coupling coefficient

        ccoefs = np.zeros((evecs.shape[1], self.nch))
        evecs = np.power(evecs, 2)

        for k in range(0, self.nch):
            ccoefs[:, k] = evecs[k*self.ngrid:(k+1)*self.ngrid, :].sum(axis=0)

        return ccoefs

    def assign_vibrational_number(self, states):

        vibns = np.empty(states.shape)
        vibns_tmp = np.empty(states.shape)
        vibns_tmp.fill(-1)

        for i, state in enumerate(states):
            vibns_tmp[state] += 1
            vibns[i] = vibns_tmp[state]

        return vibns

    def identify_levels(self):

        # select scheme for identification of calculated levels

        if self.identify == 0:
            self.get_output_identified_by_energy1()

        # slow, not recommended for large dataset
        elif self.identify == 1:
            self.get_output_identified_by_energy2()

        # fast but not recommended when exp. data contains repeated levels
        elif self.identify == 2:
            self.get_output_identified_by_coupling_coeffs1()

        elif self.identify == 3:
            self.get_output_identified_by_coupling_coeffs2()

    def view1D(self, a, b):  # a, b are arrays

        a = np.ascontiguousarray(a)
        b = np.ascontiguousarray(b)
        void_dt = np.dtype((np.void, a.dtype.itemsize * a.shape[1]))

        return a.view(void_dt).ravel(),  b.view(void_dt).ravel()

    def get_indices_of_matching_rows_argsorted(self, a, b):

        viewA, viewB = self.view1D(a, b)
        c = np.r_[viewA, viewB]
        idx = np.argsort(c, kind='mergesort')
        cs = c[idx]
        m0 = cs[:-1] == cs[1:]

        return idx[:-1][m0], idx[1:][m0]-len(viewA)

    def get_indices_of_matching_rows_searchsorted(self, a, b):

        A, B = self.view1D(a, b)
        sidxB = B.argsort()
        mask = np.isin(A, B)
        cm = A[mask]
        idx0 = np.flatnonzero(mask)
        idx1 = sidxB[np.searchsorted(B, cm, sorter=sidxB)]

        # idx0 : indices in A, idx1 : indices in B
        return idx0, idx1

    def get_indices_of_matching_rows_simple(self, a, b):

        A, B = self.view1D(a, b)
        inds = np.argwhere(A[:, None] == B)

        return inds[:, 0], inds[:, 1]

    def filter_data_by_coupling_coeffs(self):

        emask = \
            np.in1d(self.exp_data[:, 1], self.evals_predicted[:, -3]) & \
            np.in1d(self.exp_data[:, 2], self.evals_predicted[:, 1]) & \
            np.in1d(self.exp_data[:, 4], self.evals_predicted[:, 2]) & \
            np.in1d(self.exp_data[:, -1], self.evals_predicted[:, -4]) & \
            np.in1d(
                np.floor_divide(self.exp_data[:, -2], 10) + 1,
                self.evals_predicted[:, 3]
            )

        self.exp_data = self.exp_data[emask]

        cmask = \
            np.in1d(self.evals_predicted[:, -3], self.exp_data[:, 1]) & \
            np.in1d(self.evals_predicted[:, 1], self.exp_data[:, 2]) & \
            np.in1d(self.evals_predicted[:, 2], self.exp_data[:, 4]) & \
            np.in1d(self.evals_predicted[:, -4], self.exp_data[:, -1]) & \
            np.in1d(
                self.evals_predicted[:, 3],
                np.floor_divide(self.exp_data[:, -2], 10) + 1
            )

        self.evals_predicted = self.evals_predicted[cmask]

        # sort data by v then by J then by parity then by state, then by energy
        # self.exp_data = self.exp_data[
        #     np.lexsort((self.exp_data[:,3], self.exp_data[:,-1],
        #     self.exp_data[:,4], self.exp_data[:,2], self.exp_data[:,1]))
        # ]

        # self.evals_predicted = self.evals_predicted[
        #     np.lexsort((
        #             self.evals_predicted[:, 0], self.evals_predicted[:, -2],
        #             self.evals_predicted[:, 2], self.evals_predicted[:, 1],
        #             self.evals_predicted[:, -1]
        #     ))
        # ]

        i1, i2 = self.get_indices_of_matching_rows_argsorted(
            self.exp_data[:, [1, 2, 4, -1]],
            self.evals_predicted[:, [-3, 1, 2, -4]]
        )

        # self.exp_data = self.exp_data[i1]
        # self.evals_predicted = self.evals_predicted[i2]

        exp_data_ext = self.exp_data[i1]
        calc_data_ext = self.evals_predicted[i2]

        return exp_data_ext, calc_data_ext

    def filter_data_by_energy_diff(self):

        emask = \
            np.in1d(self.exp_data[:, 2], self.evals_predicted[:, 1]) & \
            np.in1d(self.exp_data[:, 4], self.evals_predicted[:, 2]) & \
            np.in1d(
                np.floor_divide(self.exp_data[:, -2], 10) + 1,
                self.evals_predicted[:, 3]
            )

        self.exp_data = self.exp_data[emask]

        cmask = \
            np.in1d(self.evals_predicted[:, 1], self.exp_data[:, 2]) & \
            np.in1d(self.evals_predicted[:, 2], self.exp_data[:, 4]) & \
            np.in1d(
                self.evals_predicted[:, 3],
                np.floor_divide(self.exp_data[:, -2], 10) + 1
            )

        self.evals_predicted = self.evals_predicted[cmask]

        # sort data by v then by J then by parity then by state, then by energy
        # self.exp_data = self.exp_data[
        #     np.lexsort((self.exp_data[:,3], self.exp_data[:,-1],
        #     self.exp_data[:,4], self.exp_data[:,2], self.exp_data[:,1]))
        # ]

        # self.evals_predicted = self.evals_predicted[
        #     np.lexsort((
        #         self.evals_predicted[:, 0], self.evals_predicted[:, -2],
        #         self.evals_predicted[:, 2], self.evals_predicted[:, 1],
        #         self.evals_predicted[:, -1]
        #     ))
        # ]

        i1, i2 = self.get_indices_of_matching_rows_simple(
            self.exp_data[:, [2, 4]],
            self.evals_predicted[: [1, 2]]
        )

        exp_data_ext = self.exp_data[i1]
        # self.evals_predicted = self.evals_predicted[i2]
        calc_data_ext = self.evals_predicted[i2]

        return exp_data_ext, calc_data_ext

    def get_output_identified_by_energy1(self):

        exp_data_ext, calc_data_ext = self.filter_data_by_energy_diff()

        coupled_coeffs = calc_data_ext[:, 4:-4]

        nlambda = np.fromiter(map(
            lambda x: self.channels[x-1].nlambda,
            np.int64(exp_data_ext[:, -1])),
            dtype=np.int64
        )

        nsigma = np.fromiter(map(
            lambda x: self.channels[x-1].nsigma,
            np.int64(exp_data_ext[:, -1])),
            dtype=np.float64
        )

        omega = np.fromiter(map(
            lambda x: self.channels[x-1].omega,
            np.int64(exp_data_ext[:, -1])),
            dtype=np.float64
        )

        out_data_ext = np.column_stack((
            exp_data_ext[:, 1],
            exp_data_ext[:, 2],
            omega, nsigma, nlambda,
            exp_data_ext[:, 4],
            exp_data_ext[:, 6],
            calc_data_ext[:, 0],
            exp_data_ext[:, 3],
            calc_data_ext[:, 0] - exp_data_ext[:, 3],
            exp_data_ext[:, 5],
            coupled_coeffs,
            exp_data_ext[:, -1],
        ))

        # sort by the absolute value of exp. and calc. energy differences
        out_data_ext = out_data_ext[np.abs(out_data_ext[:, 9]).argsort()]
        out_data_ext = out_data_ext[0:self.exp_data.shape[0], :]

        _, inds = np.unique(
            out_data_ext[:, [0, 1, 5, 8, -1]], axis=0, return_index=True
        )

        self.out_data = out_data_ext[inds]

    # Note: do not use this procedure when exp. data contains repeated data
    def get_output_identified_by_coupling_coeffs1(self):

        # self.exp_data = self.exp_data[self.exp_data[:,3].argsort()]
        # self.evals_predicted = \
        # self.evals_predicted[self.evals_predicted[:,0].argsort()]

        exp_data_ext, calc_data_ext = self.filter_data_by_coupling_coeffs()

        cstates = calc_data_ext[:, -4]
        coupled_coeffs = calc_data_ext[:, 4:-4]

        nlambda = np.fromiter(map(
            lambda x: self.channels[x-1].nlambda,
            np.int64(cstates)),
            dtype=np.int64
        )

        nsigma = np.fromiter(map(
            lambda x: self.channels[x-1].nsigma,
            np.int64(cstates)),
            dtype=np.float64
        )

        omega = np.fromiter(map(
            lambda x: self.channels[x-1].omega,
            np.int64(cstates)),
            dtype=np.float64
        )

        self.out_data = np.column_stack((
            calc_data_ext[:, -3],                       # cv
            calc_data_ext[:, 1],                        # cj
            omega, nsigma, nlambda,                     # qn
            calc_data_ext[:, 2],                        # cp
            exp_data_ext[:, 6],                         # em
            calc_data_ext[:, 0],                        # ce
            exp_data_ext[:, 3],                         # ee
            calc_data_ext[:, 0] - exp_data_ext[:, 3],   # delta
            exp_data_ext[:, 5],                         # eu
            coupled_coeffs,                             # cc
            cstates                                     # cst
        ))

    def get_output_identified_by_coupling_coeffs2(self):

        # self.exp_data = self.exp_data[self.exp_data[:, 3].argsort()]
        # self.evals_predicted = \
        # self.evals_predicted[self.evals_predicted[:, 0].argsort()]

        self.filter_data_by_coupling_coeffs()

        used_data = []

        for expd in self.exp_data[:, 1:]:
            ev, ej, ee, ep, eu, em, es = expd

            for cn, calcd in enumerate(self.evals_predicted[:, :-2]):

                ce, cj, cp, iso = calcd[0:4]
                cc = calcd[4:-2]
                cs, cv = calcd[-2:]

                if cj == ej and cp == ep and np.int64(cv) == np.int64(ev) and \
                   np.int64(cs) == np.int64(es) and \
                   np.int64(iso) == np.int64(em / 10) + 1 and \
                   cn not in used_data:

                    vib, rot = np.int64(cv), cj
                    args = (np.int64(cp), np.int64(em), ce, ee, ce-ee, eu)
                    st = np.int64(cs)
                    lam = np.int64(self.channels[st-1].nlambda)
                    sig = np.float64(self.channels[st-1].nsigma)
                    omg = np.float64(self.channels[st-1].omega)

                    line = np.array([vib, rot, omg, sig, lam, *args, *cc, st])

                    used_data.append(cn)

                    self.out_data = np.append(self.out_data, line)
                    break

        self.out_data = self.out_data.reshape(-1, 12+self.nch)

    # Note: very slow
    def get_output_identified_by_energy2(self):

        # self.exp_data = self.exp_data[self.exp_data[:, 3].argsort()]
        # self.evals_predicted = \
        # self.evals_predicted[self.evals_predicted[:, 0].argsort()]
        # start1 = time.time()

        # self.filter_data_for_energy_identification()

        used_data = []

        for expd in self.exp_data[:, 1:]:

            ev, ej, ee, ep, eu, em, es = expd

            vib, rot, st = -1, -1.0, -1
            lam, sig, omg = -1, -1.0, -1.0
            args, ccoefs = None, None
            min_diffe = 1.0e10
            is_diffe_min = False
            cnum = 0

            for calcn, calcd in enumerate(self.evals_predicted[:, :-2]):
                ce, cj, cp, iso = calcd[0:4]
                cc = calcd[4:-2]

                diffe = abs(ee - ce)

                if cj == ej and cp == ep and \
                   np.int64(iso) == np.int64(em / 10) + 1:
                    # and ev <= self.emaxi:  #and calcn not in used_data:

                    if diffe < min_diffe:
                        is_diffe_min = False
                        min_diffe = diffe

                        vib, rot = np.int64(ev), cj
                        args = (np.int64(cp), np.int64(em), ce, ee, ce-ee, eu)
                        st = np.int64(es)
                        lam = np.int64(self.channels[st-1].nlambda)
                        sig = np.float64(self.channels[st-1].nsigma)
                        omg = np.float64(self.channels[st-1].omega)
                        ccoefs = cc
                    else:
                        is_diffe_min = True
                        cnum = calcn
                        break

            if is_diffe_min:
                # line will be an array of type object not np.float64
                used_data.append(cnum)
                line = np.array([vib, rot, omg, sig, lam, *args, *ccoefs, st])

                self.out_data = np.append(self.out_data, line)

        self.out_data = self.out_data.reshape(-1, 12+self.nch)

    def _save_eigenvectors(self, j, par, iso, evector, vibnums, states):

        # TODO: write as binary files

        os.makedirs(self.evec_dir, exist_ok=True)
        fname = f'evector_J{str(math.floor(j))}_p{str(par)}_i{str(iso)}.dat'
        efile = os.path.join(self.evec_dir, fname)

        evec_extended = np.vstack((vibnums, states, evector))

        # np.savetxt(efile,
        # np.column_stack((self.rgrid, evector)), newline='\n', fmt='%15.8e')
        np.savetxt(efile, evec_extended, newline='\n', fmt='%15.8e')

    # TODO: Is this needed here?
    def interpolate_wavefunction(self, ninter=200, jrange=(0, 1), pars=(1,),
                                 vrange=(0, 1), iso=(1,)):

        igrid, istep = np.linspace(
            self.rmin,
            self.rmax,
            num=ninter,
            endpoint=False,  # ??
            retstep=True
        )

        efiles = [
            f'evector_J{j}_{k}_{i}.dat'
            for j in jrange for k in pars for i in iso
        ]

        for efile in efiles:  # os.listdir(self.evec_dir):
            nvib = np.arange(vrange[0], vrange[1] + 1, dtype=np.int64)
            evec = np.loadtxt(
                os.path.join(self.evec_dir, efile),
                usecols=tuple(nvib+1)
            )  # skip zero column
            arr = []

            for col in nvib:  # range(0, nvib):
                for row in range(0, ninter):
                    index = 0
                    sum_arg = 0.0
                    for j in range(0, self.nch*self.ngrid):

                        if j % self.ngrid == 0:
                            index += 1

                        expr = igrid[row] - self.rgrid[j-index*self.ngrid]
                        if abs(expr) < 1e-15:
                            sinc = 1.0
                        else:
                            arg = (math.pi / istep) * expr
                            sinc = math.sin(arg) / arg

                        sum_arg += evec[j][col] * sinc

                    arr.append(sum_arg)

            wavefunc = np.array(arr).reshape(len(nvib), ninter).T
            nwavefunc = wavefunc / wavefunc.sum(axis=0, keepdims=1)

            self.save_wavefunctions(np.column_stack((igrid, nwavefunc)))

    def save_wavefunctions(self, wavefunc):

        os.makedirs(self.wavef_dir, exist_ok=True)

        efile = os.path.join(self.wavef_dir, 'wavefunctions.dat')

        fmt = wavefunc.shape[1] * ['%15.8e']

        np.savetxt(efile, wavefunc, newline='\n', fmt=fmt)

    def calculate_stats(self, yobs, ycal, yunc, is_weighted):

        diff_square = np.square(yobs - ycal)

        # calculate chi square
        if not is_weighted:
            weights = 1.0 / np.square(yunc)
        else:
            weights = 1.0 / (np.square(yunc) + 0.33 * (diff_square))

        chi2 = np.sum(diff_square * weights) / yobs.shape[0]

        # calculate rms
        s = np.sum(diff_square) / yobs.shape[0]
        rms = math.sqrt(s)

        # calculate dimensionless rms
        rmsd = math.sqrt(chi2)

        return chi2, rms, rmsd

    def save_output_data(self, append, stats):

        cc = 'CC'
        cc_names = ''.join(map(
            lambda x: f'{cc+str(x):>9}', range(1, self.nch+1))
        )

        names = [
            'No', 'v', 'J', 'omega', 'sigma',
            'lambda', 'parity', 'marker',
            'Ecalc', 'Eexp', 'delta',
            'unc', cc_names, 'state'
        ]

        header = '{:>6}'.format(names[0]) + \
            ''.join(map(lambda x: '{:>8}'.format(x), names[1:9])) + \
            ''.join(map(lambda x: '{:>15}'.format(x), names[9:11])) + \
            ''.join(map(lambda x: '{:>11}'.format(x), names[11:]))

        self.out_data = np.c_[
            np.arange(self.ninit, self.out_data.shape[0]+self.ninit),
            self.out_data
        ]

        fmt = [
            '%7d', '%7d', '%7.1f', '%7.1f', '%7.1f', '%6d', '%6d',
            '%6d', '%15.6f', '%14.6f', '%12.6f', '%10.4f'
        ] + self.nch*['%8.3f'] + ['%5d']

        footer = \
            f'Chi Square = {stats[0]:.8f}\n' \
            f'RMS = {stats[1]:.8f} cm-1\n' \
            f'RMSD = {stats[2]:.8f}'

        wmode = 'w'
        if append:
            wmode = 'a'
            header = footer = ''

        if len(self.masses) > 1:
            footer = ''

        with open(self.file_evals, wmode) as outs:
            np.savetxt(
                outs,
                self.out_data,
                comments='#',
                header=header,
                footer=footer,
                fmt=fmt
            )

    def save_all_calculated_data(self):

        cols = [-3] + list(
            range(0, self.evals_predicted.shape[1]-3)
        ) + [-2, -1]

        eigv_data = self.evals_predicted[:, cols]
        eigv_data = np.c_[
            np.arange(1.0, self.evals_predicted.shape[0]+1), eigv_data
        ]

        cc = 'CC'
        cc_names = ''.join(map(
            lambda x: f'{cc+str(x):>9}', range(1, self.nch+1))
        )

        names = [
            'No', 'v', 'Ecalc', 'J', 'parity', ' marker',
            cc_names, 'state', 'lambda', 'omega'
        ]

        header = \
            '{:>7}'.format(names[0]) + \
            ''.join(map(lambda x: '{:>4}'.format(x), names[1])) + \
            ''.join(map(lambda x: '{:>11}'.format(x), names[2:4])) + \
            ''.join(map(lambda x: '{:>8}'.format(x), names[4:]))

        fm = ['%7.1d', '%3.1d', '%14.6f', '%7.1f', '%5.1d', '%6.1d'] + \
            self.nch*['%9.3f'] + 2*['%6.1d'] + ['%8.2f']

        np.savetxt(
            self.file_predicted, eigv_data, comments='#', header=header, fmt=fm
        )

    def save_additional_data(self):

        sp = 4*' '

        with open(self.info_outfile, 'w') as outf:
            outf.write(f'{30*"#"} Grid and Molecule Data {30*"#"}\n\n')
            outf.write(f'{sp}Number of Grid Points = {self.ngrid:<5d}\n\n')
            outf.write(
                f'{sp}Rmin = {self.rmin*Const.bohr:>10.8f} '
                f'Angstrom = {self.rmin:>10.8f} Bohr\n'
            )
            outf.write(
                f'{sp}Rmax = {self.rmax*Const.bohr:>10.8f} '
                f'Angstrom = {self.rmax:>10.8f} Bohr\n\n'
            )

            rgrid_str = np.array2string(self.rgrid)
            outf.write(f'{sp}Grid Points:\n{sp}{rgrid_str}\n\n')

            outf.write(f'{sp}Method of solution: {self.solver}\n')

            outf.write(f'\n{sp}Number of Isotopes = {len(self.masses):>15d}\n')
            for im, mass in enumerate(self.masses):
                outf.write(
                    f' {sp}{im+1}. {mass:>10.9} au = '
                    f'{mass/Const.massau:>10.9} amu\n'
                )

            outf.write(f'\n{sp}Rotational Qunatum Numbers = ')
            outf.write(', '.join(map(str, self.jqnumbers)))
            outf.write(
                f'\n{sp}Number of Calculated Eigenvalues = '
                f'{len(self.evals_predicted):>15d}'
            )
            outf.write(
                f'\n{sp}Number of Selected Eigenvalues = '
                f'{len(self.out_data):>15d}'
            )
            # outf.write(f'\n{sp}Maximum Number of Eigenvalues = {self.emaxi}')

            ofj = True if self.refj is not None else False
            outf.write(f'\n{sp}Shift Energies = {ofj}')
            if ofj:
                outf.write(f' by Level J = {self.refj}')

            outf.write(f'\n{sp}Levels with Parity = ')
            outf.write(', '.join(map(lambda x: 'e' if x else 'f', self.pars)))

            outf.write(
                f'\n{sp}File with Experimental Data --- {self.exp_file}\n'
            )

            outf.write(f'\n{30*"#"} Channels Data {30*"#"}\n')
            outf.write(f'\n{sp}Number of Channels = {self.nch}\n')

            for ic, ch in enumerate(self.channels):
                outf.write(f'\n{sp:>5} {ic+1}. Type = {ch.model:>10}\n')
                outf.write(f'{sp:>9}Lambda = {ch.nlambda:>3}\n')
                outf.write(f'{sp:>9}Sigma = {ch.nsigma:>3}\n')
                outf.write(f'{sp:>9}Multiplicity = {ch.mult:>3}\n')
                outf.write(f'{sp:>9}Parameters:\n')
                if ch.model == 'pointwise':
                    for i in range(0, len(ch.R)):
                        outf.writelines(
                            f'{sp} {ch.R[i]:>15.10f}{sp}'
                            f'{ch.U[i]*Const.hartree:>15.10f}\n'
                        )

            outf.write(f'\n{30*"#"} Couplings Data {30*"#"}\n')
            outf.write(f'\n{sp}Number of Couplings = {self.ncp}')
            outf.write(
                f'\n{sp}Hamiltonian Matrix Size = '
                f'{self.nch*self.ngrid}x{self.nch*self.ngrid}\n'
            )
            outf.write(f'\n{30*"#"} Coupling Functions on Grid {30*"#"}\n\n')

        fgrid_cols = np.hstack((
            self.rgrid[:, np.newaxis],
            self.fgrid.reshape(self.ncp, self.ngrid).T
        ))

        # enames = 'N', 'v', 'Ecalc', 'J', 'parity', 'marker', 'state'
        # eheader = '{:>7}'.format(enames[0]) + \
        #     ''.join(map(lambda x: '{:>6}'.format(x), enames[1])) + \
        #     ''.join(map(lambda x: '{:>14}'.format(x), enames[2:3])) + \
        #     ''.join(map(lambda x: '{:>11}'.format(x), enames[3:4])) + \
        #     ''.join(map(lambda x: '{:>8}'.format(x), enames[4:]))

        # efmt = ['%7.1d', '%3.1d', '%14.6f', '%7.1f', '%4.1d', '%4.1d']

        with open(self.info_outfile, 'a') as outf:
            np.savetxt(outf, fgrid_cols, fmt=fgrid_cols.shape[1] * ['%16.10f'])
            outf.write(f'\n\n{20*"#"} Experimental data {20*"#"}\n\n')
            outf.write(
                f'\nNumber of used experimental data = {len(self.exp_data)}\n'
            )
            # np.savetxt(
            #   outf, self.exp_data, comments='', header=eheader, fmt=efmt
            # )
            # outf.write(f'\n\n{20*"#"} Eigenvalues {20*"#"}\n\n')
            # np.savetxt(outf, eigv_data, comments='', header=header, fmt=fmt)
