import math
import os
import numpy as np
import scipy as sp
# from numba import njit
from scipy.interpolate import CubicSpline
from .interaction import Interaction
from .molecule_data import Channel, Coupling
from constants import Const
from .grids import CSpline
from .utils import Utils
import Utils.C_hartree as C_hartree
import Utils.C_bohr as C_bohr
import Utils.C_massau as C_massau


class MoleculeLevels:

    def __init__(self, md, grid, channels, couplings=[], data=None):

        self.jqnumbers = md.jqnumbers
        self.pars = md.parities
        self.masses = md.reduced_masses or md.masses
        self.nisotopes = md.nisotopes
        self.refj = md.referencej
        self.refE = md.refE
        self.exp_data = md.exp_data

        # filter by state
        state_numbers = np.arange(1, len(channels)+1, dtype=np.int64)
        state_mask = np.in1d(
            self.exp_data[:, -1], np.fromiter(state_numbers, dtype=np.float64)
        )
        self.exp_data = self.exp_data[state_mask]

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

        self.eig_decomp_keys = self._get_eig_decomp_keys()

        self.solver_keys = self._get_solver_keys()

        self.rgrid2 = np.square(self.rgrid)
        self.ugrid = np.zeros(self.nch * self.ngrid)
        self.fgrid = np.zeros(self.ncp * self.ngrid)

        # spline derivatives
        self.sderiv = False

        # maximum number of fit parameters
        nmax_params = 200
        self.sk_grid = np.zeros((self.nch * self.ngrid, nmax_params))

        self.pot_enr = np.zeros((self.nch*self.ngrid, self.nch*self.ngrid))
        self.kin_enr = np.zeros_like(self.pot_enr)

        # initilizie the total number of potential parameters
        self._ctot = 0

        self.evals_file = 'eigenvalues.dat'
        self.predicted_file = 'evalues_predicted.dat'
        self.info_outfile = 'data_info.dat'

        self.evec_dir = os.path.join(Utils.get_current_dir(), 'eigenvectors')
        self.wavef_dir = os.path.join(Utils.get_current_dir(), 'wavefunctions')

    def _get_eig_decomp_keys(self):

        return {
            'lapack': self._lapack_eig_decomposition,
            'arpack': self._arpack_eig_decomposition
        }

    def _get_solver_keys(self):

        return {
            'sinc': self._calculate_kinetic_energy_sinc,
            'fourier': self._calculate_kinetic_energy_fourier,
            'fd5': self._calculate_kinetic_energy_FD5
        }

    def calculate_channels_on_grid(self, ypar=None):

        self._ctot = 0
        unique = set()
        pot_ranges = dict()

        # pass by reference
        index = [-1]

        for ch in range(1, self.nch+1):

            # determine _ctot value
            if self.channels[ch-1].filep not in unique:
                unique.add(self.channels[ch-1].filep)

                self._ctot += self.channels[ch-1].npnts

            # built-in cubic spline function
            if self.channels[ch-1].model.lower() == Channel.models[1]:
                self._calculate_pointwise_pec_on_grid(
                    ch, ypar, pot_ranges, index
                )

            # cubic spline function using with custom implementation
            if self.channels[ch-1].model.lower() == Channel.models[2]:
                self._calculate_custom_pointwise_pec_on_grid(
                    ch, ypar, pot_ranges, index
                )

            # Morse potential function
            if self.channels[ch-1].model.lower() == Channel.models[3]:
                self._calculate_Morse_pec_on_grid(ch, ypar, pot_ranges, index)

            # EMO potential function
            if self.channels[ch-1].model.upper() == Channel.models[4]:
                self._calculate_EMO_pec_on_grid(ch, ypar, pot_ranges, index)

            # MLR potential function
            if self.channels[ch-1].model.upper() == Channel.models[5]:
                self._calculate_MLR_pec_on_grid(ch, ypar, pot_ranges, index)

            # Custom potential function
            if self.channels[ch-1].model.lower() == Channel.models[6]:
                self._calculate_custom_pec_on_grid(ch, ypar, pot_ranges, index)

    def _calculate_pointwise_pec_on_grid(self, ch, ypar, pot_ranges, index):

        # xpnts, ypnts = self.calculate_pec_points(ch, ypar, pot_ranges, index)

        xpnts = self.channels[ch-1].R

        # ypar contains parameters from all potentials and they
        # have to be transformed to each channel
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

    def _calculate_custom_pointwise_pec_on_grid(self, ch, ypar, pot_ranges,
                                                index):

        xpnts = self.channels[ch-1].R
        ypnts, st, en = self._calculate_pec_points(ch, ypar, pot_ranges, index)

        cs = CSpline(xpnts, ypnts)
        index = (ch-1)*self.ngrid, ch*self.ngrid

        self.ugrid[index[0]:index[1]], sk = \
            cs(self.rgrid, return_deriv=True)

        self.sk_grid[index[0]:index[1], st:en] = sk
        self.sderiv = True

    def _calculate_pec_points(self, ch, ypar, pot_ranges, index):

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
            st = (ch-1)*self.channels[ch-1].npnts
            en = ch*self.channels[ch-1].npnts
            ypnts = self.channels[ch-1].U

        return ypnts, st, en

    def _calculate_Morse_pec_on_grid(self, ch, ypar, pot_ranges, index):

        ypnts, *_ = self._calculate_pec_points(ch, ypar, pot_ranges, index)

        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
            self._VMorse(ypnts, self.rgrid)

    def _VMorse(self, params, rgrid):

        Te, De, a, re = params

        return Te + De*np.power((1.0 - np.exp(-a*(rgrid-re))), 2.0)

    def _calculate_EMO_pec_on_grid(self, ch, ypar, pot_ranges, index):

        if ypar is not None:
            pass
        else:
            self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
                self._VEMO(self.channels[ch-1].U, self.rgrid)

    def _VEMO(self, params, rgrid):

        Te, De, p, re = params[:4]
        bparams = np.array(params[4:])[::-1]

        yr = (np.power(rgrid, p) - re ** p) / (np.power(rgrid, p) + re ** p)

        bemo = np.polyval(bparams, yr)
        vemo = Te + De * np.power((1.0 - np.exp((-1.0*bemo)*(rgrid-re))), 2.0)

        return vemo

    def _calculate_MLR_pec_on_grid(self, ch, ypar, pot_ranges, index):

        if ypar is not None:
            pass
        else:
            ni = self.channels[ch-1].ni
            nb = self.channels[ch-1].nb
            nc = self.channels[ch-1].nc
            nd = self.channels[ch-1].nd

            self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
                self._VMLR(self.channels[ch-1].U, self.rgrid, (ni, nb, nc, nd))

    def _VMLR(self, params, rgrid, nparams):

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

    def _calculate_custom_pec_on_grid(self, ch, ypar, pot_ranges, index):

        ypnts, *_ = self._calculate_pec_points(ch, ypar, pot_ranges, index)

        self.ugrid[(ch-1)*self.ngrid:ch*self.ngrid] = \
            self.channels[ch-1].cfunc(ypnts, self.rgrid)

    def get_potential_on_grid(self):

        return self.ugrid

    def get_grid_points(self):

        return self.rgrid

    def get_couplings_on_grid(self):

        return self.fgrid

    def calculate_couplings_on_grid(self, ypar=None):

        # pass by reference
        yrange = {'start': self._ctot, 'end': self._ctot}

        for cp in range(0, self.ncp):

            # built-in cubic spline function
            if self.couplings[cp].model.lower() == Coupling.models[1]:
                self._calculate_pointwise_couplings_on_grid(cp, yrange, ypar)

            # custom implementation of cubic spline function
            if self.couplings[cp].model.lower() == Coupling.models[2]:
                self._calculate_custom_pointwise_couplings_on_grid(
                    cp, yrange, ypar
                )

            # constant coupling value
            if self.couplings[cp].model.lower() == Coupling.models[3]:
                self.calculate_constant_coupling_on_grid(cp, ypar)

            # custom coupling function
            if self.couplings[cp].model.lower() == Coupling.models[4]:
                self._calculate_custom_coupling_on_grid(cp, yrange, ypar)

    def _calculate_pointwise_couplings_on_grid(self, cp, yrange, ypar):

        xpnts = self.couplings[cp].xc
        ypnts = self.calculate_coupling_points(cp, yrange, ypar)

        cs = CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)

        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = cs(self.rgrid)

    def _calculate_custom_pointwise_couplings_on_grid(self, cp, yrange, ypar):

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

    def _calculate_custom_coupling_on_grid(self, cp, yrange, ypar):

        ypnts = self.calculate_coupling_points(cp, yrange, ypar)

        self.fgrid[cp*self.ngrid:(cp+1)*self.ngrid] = \
            self.couplings[cp].cfunc(ypnts, self.rgrid)

    def _calculate_kinetic_energy(self, mass):

        self.solver_keys[self.solver](mass)

    def _lapack_eig_decomposition(self):

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

    def _arpack_eig_decomposition(self):

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

    def _get_levels_keys(self):

        return {
            1: 'eig_decomp',
            2: 'energy_subset_index',
            3: 'energy_subset_value',
            4: 'identify',
            5: 'evecs_binary',
            7: 'store_evecs',
            8: 'weighted',
            11: 'lapack_driver',
            12: 'arpack_which',
            13: 'arpack_k',
            14: 'arpack_sigma',
            15: 'file_evals_predicted',
            16: 'file_eigenvalues',
        }

    def calculate_levels(self, ypar=None, **kwargs):

        # ypar is only for internal use

        keys = self._get_levels_keys()

        self.eig_decomp = kwargs.get(keys[1]) or 'lapack'
        self.energy_subset_index = kwargs.get(keys[2])
        self.energy_subset_value = kwargs.get(keys[3])
        self.identify = kwargs.get(keys[4]) or 0
        self.evecs_binary = kwargs.get(keys[5]) or True
        self.store_evecs = kwargs.get(keys[7]) or False
        self.is_weighted = kwargs.get(keys[8]) or False
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
            self.eminv = self.energy_subset_value[0] / C_hartree
            self.emaxv = self.energy_subset_value[-1] / C_hartree
            self.evalues_subsetv = (self.eminv, self.emaxv)

        chi2, rms, _ = self.calculate_eigenvalues(ypar)

        print(f'Chi Square = {chi2:<18.8f} | RMS = {rms:<15.8f}cm-1')
        # f'\nRMSD = {rmsd:.8f}'

    def calculate_eigenvalues(self, ypar=None):

        self.calculate_channels_on_grid(ypar=ypar)

        if self.ncp != 0:
            self.calculate_couplings_on_grid(ypar=ypar)

        interact = Interaction(self.ngrid, self.nch, self.rgrid2)
        dd = np.diag_indices(self.ngrid)
        count = 0

        for niso, iso in enumerate(self.nisotopes):

            self.evals_predicted = np.zeros(
                (2 * self.nch * self.exp_data.shape[0], 6+self.nch)
            )

            mass = self.masses[iso-1]

            self._calculate_kinetic_energy(mass)

            for par in self.pars:
                # pass by reference
                shiftE = [0.0]

                for jrotn in self.jqnumbers:

                    self._calculate_potential_energy(jrotn, par, mass, dd)

                    self.hmatrix = self.kin_enr + self.pot_enr

                    pert_matrix = interact.get_perturbations_matrix(
                        jrotn, mass, par, self.channels, self.couplings,
                        self.fgrid, dd, self.Gy2
                    )
                    # pmatrix = pmatrix*np.tile(self.Gy2,(pmatrix.shape[0],5))

                    self.hmatrix += pert_matrix

                    evalues, evecs = self.eig_decomp_keys[self.eig_decomp]()

                    ccoefs = self._get_coupling_coefficients(evecs)
                    states = np.argmax(ccoefs, axis=1) + 1
                    vibnums = self._assign_vibrational_number(states)

                    if self.store_evecs:

                        self._save_eigenvectors(
                            jrotn, par, iso, evecs[:, 0:evalues.shape[0]],
                            vibnums, states
                        )

                    calc_data = self._arange_levels(
                        jrotn, par, iso, shiftE, evalues
                    )

                    calc_data = np.column_stack(
                        (calc_data, ccoefs, states, vibnums)
                    ).reshape(evalues.shape[0], -1)

                    nrow, ncol = calc_data.shape[0], calc_data.shape[1]
                    self.evals_predicted[count:count+nrow, 0:ncol] = calc_data

                    count += calc_data.shape[0]

            nlambda = self._get_lambda_values()
            omega = self._get_omega_values()

            self.evals_predicted = np.column_stack(
                (self.evals_predicted, nlambda, omega)
            )

            if self.exp_data is not None:

                if niso == 0:
                    self.out_data = self._identify_levels()
                else:
                    self.out_data = np.vstack(
                        (self.out_data, self._identify_levels())
                    )
            else:
                self.save_predicted()

        chi2, rms, rmsd = self.calculate_stats(
            self.out_data[:, 8], self.out_data[:, 7],
            self.out_data[:, 10], self.is_weighted
        )

        self._save_output_data(stats=(chi2, rms, rmsd))

        return chi2, rms, rmsd

    def _get_lambda_values(self):

        return np.fromiter(map(
            lambda x: self.channels[x-1].nlambda,
            np.int64(self.evals_predicted[:, -2])),
            dtype=np.int64
        )

    def _get_omega_values(self):

        return np.fromiter(map(
            lambda x: self.channels[x-1].omega,
            np.int64(self.evals_predicted[:, -2])),
            dtype=np.float64
        )

    def _arange_levels(self, jrotn, par, iso, shiftE, evalues):

        if jrotn == self.jqnumbers[0] and self.refj is not None:
            shiftE[0] = evalues[0]

        elif self.refE is not None:
            shiftE[0] = self.refE

        evalues = (evalues - shiftE[0]) * C_hartree

        calc_data = np.column_stack((
                evalues[:evalues.shape[0]],
                np.full((evalues.shape[0],), jrotn),
                np.full((evalues.shape[0],), par),
                np.full((evalues.shape[0],), iso),
        ))

        return calc_data

    def _calculate_kinetic_energy_fourier(self, mass):

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

    def _calculate_kinetic_energy_sinc(self, mass):

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

    def _calculate_kinetic_energy_FD5(self, mass):

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
            (h2/(step**2)) * self._form_five_diagonal_symm_matrix(d0, d1, d2)

        corner_coef = 29.0 / 12.0
        self.kin_enr[0, 0] = corner_coef
        self.kin_enr[self.ngrid-1, self.ngrid-1] = corner_coef

        for ch in range(2, self.nch+1):
            ind3 = (ch-1)*self.ngrid
            ind4 = ch*self.ngrid

            self.kin_enr[ind3:ind4, ind3:ind4] = \
                self.kin_enr[0:self.ngrid, 0:self.ngrid]

    def _form_five_diagonal_symm_matrix(self, a, b, c, k1=0, k2=1, k3=2):

        return np.diag(a, k1) + np.diag(b, k2) + np.diag(b, -k2) + \
                np.diag(c, k3) + np.diag(c, -k3)

    def _calculate_potential_energy(self, jrotn, par, mass, dd):

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

            self.pot_enr[i1:i2, i1:i2][dd] = \
                self.Gy2 * (self.ugrid[i1:i2] + (num / denom))

    def _get_coupling_coefficients(self, evecs):
        # assign state based on the largest coupling coefficient

        ccoefs = np.zeros((evecs.shape[1], self.nch))
        evecs = np.power(evecs, 2)

        for k in range(0, self.nch):
            ccoefs[:, k] = evecs[k*self.ngrid:(k+1)*self.ngrid, :].sum(axis=0)

        return ccoefs

    def _assign_vibrational_number(self, states):

        vibns = np.empty(states.shape[0])
        vibns_tmp = np.empty(states.shape[0]+1)
        vibns_tmp.fill(-1)

        for i, state in enumerate(states):
            vibns_tmp[state] += 1
            vibns[i] = vibns_tmp[state]

        return vibns

    def _identify_levels(self):

        # fast identification by energy
        if self.identify == 0:
            return self._get_output_identified_by_energy1()

        # fast identification by coupling cofficients
        elif self.identify == 1:
            return self._get_output_identified_by_coupling_coeffs1()

        # relatively fast identification by coupling cofficients
        elif self.identify == 2:
            return self._get_output_identified_by_coupling_coeffs2()

        # very slow identification by energy, used for debugging
        elif self.identify == 3:
            return self._get_output_identified_by_energy2()

    @staticmethod
    def _view1D(a, b):  # a, b are arrays

        a = np.ascontiguousarray(a)
        b = np.ascontiguousarray(b)
        void_dt = np.dtype((np.void, a.dtype.itemsize * a.shape[1]))

        return a.view(void_dt).ravel(),  b.view(void_dt).ravel()

    @staticmethod
    def _get_indices_of_matching_rows_argsorted(a, b):

        viewA, viewB = MoleculeLevels._view1D(a, b)
        c = np.r_[viewA, viewB]
        idx = np.argsort(c, kind='mergesort')
        cs = c[idx]
        m0 = cs[:-1] == cs[1:]

        return idx[:-1][m0], idx[1:][m0]-len(viewA)

    @staticmethod
    def _get_indices_of_matching_rows_searchsorted(a, b):

        A, B = MoleculeLevels._view1D(a, b)
        sidxB = B.argsort()
        mask = np.isin(A, B)
        cm = A[mask]
        idx0 = np.flatnonzero(mask)
        idx1 = sidxB[np.searchsorted(B, cm, sorter=sidxB)]

        # idx0 : indices in A, idx1 : indices in B
        return idx0, idx1

    @staticmethod
    def _get_indices_of_matching_rows_simple(a, b):

        A, B = MoleculeLevels._view1D(a, b)
        inds = np.argwhere(A[:, None] == B)

        return inds[:, 0], inds[:, 1]

    def _filter_data_by_coupling_coeffs(self):

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
        # self.evals_predicted = self.evals_predicted[np.lexsort((
        #             self.evals_predicted[:, 0], self.evals_predicted[:, -2],
        #             self.evals_predicted[:, 2], self.evals_predicted[:, 1],
        #             self.evals_predicted[:, -1]
        # ))]

        i1, i2 = self._get_indices_of_matching_rows_simple(
            self.exp_data[:, [1, 2, 4, -1]],
            self.evals_predicted[:, [-3, 1, 2, -4]]
        )

        exp_data_ext = self.exp_data[i1]
        calc_data_ext = self.evals_predicted[i2]

        return exp_data_ext, calc_data_ext

    def _filter_data_by_energy_diff(self):

        emask = \
            np.in1d(self.exp_data[:, 2], self.evals_predicted[:, 1]) & \
            np.in1d(self.exp_data[:, 4], self.evals_predicted[:, 2]) & \
            np.in1d(
                np.floor_divide(self.exp_data[:, -2], 10) + 1,
                self.evals_predicted[:, 3]
            )

        exp_data_mask = self.exp_data[emask]

        cmask = \
            np.in1d(self.evals_predicted[:, 1], exp_data_mask[:, 2]) & \
            np.in1d(self.evals_predicted[:, 2], exp_data_mask[:, 4]) & \
            np.in1d(
                self.evals_predicted[:, 3],
                np.floor_divide(exp_data_mask[:, -2], 10) + 1
            )

        self.evals_predicted = self.evals_predicted[cmask]

        # sort data by v then by J then by parity then by state, then by energy
        # self.exp_data = self.exp_data[
        #     np.lexsort((self.exp_data[:,3], self.exp_data[:,-1],
        #     self.exp_data[:,4], self.exp_data[:,2], self.exp_data[:,1]))
        # ]
        # self.evals_predicted = self.evals_predicted[np.lexsort((
        #         self.evals_predicted[:, 0], self.evals_predicted[:, -2],
        #         self.evals_predicted[:, 2], self.evals_predicted[:, 1],
        #         self.evals_predicted[:, -1]
        # ))]
        i1, i2 = self._get_indices_of_matching_rows_simple(
            exp_data_mask[:, [2, 4]],
            self.evals_predicted[:, [1, 2]]
        )

        exp_data_ext = exp_data_mask[i1]
        calc_data_ext = self.evals_predicted[i2]

        return exp_data_ext, calc_data_ext

    def _get_output_identified_by_energy1(self):

        exp_data_ext, calc_data_ext = self._filter_data_by_energy_diff()
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

        _, inds = np.unique(
            out_data_ext[:, [0, 1, 5, 8, -1]], axis=0, return_index=True
        )

        return out_data_ext[inds]

    def _get_output_identified_by_coupling_coeffs1(self):

        # self.exp_data = self.exp_data[self.exp_data[:,3].argsort()]
        # self.evals_predicted = \
        # self.evals_predicted[self.evals_predicted[:,0].argsort()]

        exp_data_ext, calc_data_ext = self._filter_data_by_coupling_coeffs()

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

        out_data_ext = np.column_stack((
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

        return out_data_ext

    def _get_output_identified_by_coupling_coeffs2(self):

        # self.exp_data = self.exp_data[self.exp_data[:, 3].argsort()]
        # self.evals_predicted = \
        # self.evals_predicted[self.evals_predicted[:, 0].argsort()]

        self._filter_data_by_coupling_coeffs()

        used_data = []
        out_data_ext = np.array([], dtype=np.float64)

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

                    out_data_ext = np.append(out_data_ext, line)
                    break

        return out_data_ext.reshape(-1, 12+self.nch)

    def _get_output_identified_by_energy2(self):

        # self.exp_data = self.exp_data[self.exp_data[:, 3].argsort()]
        # self.evals_predicted = \
        # self.evals_predicted[self.evals_predicted[:, 0].argsort()]
        # self.filter_data_for_energy_identification()

        used_data = []
        out_data_ext = np.array([], dtype=np.float64)

        for expd in self.exp_data[:, 1:]:

            ev, ej, ee, ep, eu, em, es = expd
            vib, rot, st = -1, -1.0, -1
            lam, sig, omg = -1, -1.0, -1.0
            min_diffe = 1.0e10
            cnum = 0
            args, ccoefs = None, None
            is_diffe_min = False

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

                out_data_ext = np.append(out_data_ext, line)

        return out_data_ext.reshape(-1, 12+self.nch)

    def _save_eigenvectors(self, jn, par, iso, evector, vibnums, states):

        os.makedirs(self.evec_dir, exist_ok=True)
        fname = f'evector_J{str(math.floor(jn))}_p{str(par)}_i{str(iso)}'
        efile = os.path.join(self.evec_dir, fname)

        evec_ext = np.vstack((vibnums, states, evector))

        if self.evecs_binary:
            np.save(efile, evec_ext)
        else:
            efile = efile + '.dat'
            np.savetxt(efile, evec_ext, newline='\n', fmt='%19.12e')

    def interpolate_wavefunction(self, ninter=200, jrange=(0, 1),
                                 pars=(1,), binary=True):
        """Interpolate the stored eigenvectors using sinc basis functions

        Args:
            ninter (int, optional): Number of interpolation points.
            jrange (tuple, optional): The range of rotational quantum numbers.
            pars (tuple, optional): The parities which to be used.
            binary (bool, optional): store the wavefunction as binary file.
        """

        igrid, istep = np.linspace(
            self.rmin,
            self.rmax,
            num=ninter,
            endpoint=False,
            retstep=True
        )

        efiles = [
            os.path.join(self.evec_dir, f'evector_J{j}_p{k}_i{i}.dat')
            for j in jrange for k in pars for i in self.nisotopes
        ]

        gridr = np.tile(self.rgrid, self.nch)

        for efile in efiles:
            evec = np.loadtxt(efile)
            for col in range(evec.shape[1]):
                wavefunc = self.sinc_interp(evec[2:, col], gridr, igrid)
                nwavefunc = wavefunc / np.linalg.norm(wavefunc)

                self.save_wavefunctions(np.column_stack((igrid, nwavefunc)))

    def sinc_interp(self, x, s, u):
        """
        Interpolates x, sampled at "s" instants
        Output y is sampled at "u" instants ("u" for "upsampled")
        """

        if len(x) != len(s):
            raise ValueError('x and s must be the same length')

        # Find the period
        T = s[1] - s[0]

        sinc_m = np.tile(u, (len(s), 1))-np.tile(s[:, np.newaxis], (1, len(u)))
        return np.dot(x, np.sinc(sinc_m/T))

    # def interpolate_wavefunction(self, ninter=200, jrange=(0, 1), pars=(1,)):
    #     igrid, istep = np.linspace(
    #         self.rmin,
    #         self.rmax,
    #         num=ninter,
    #         endpoint=False,  # ??
    #         retstep=True
    #     )
    #     efiles = [
    #         os.path.join(self.evec_dir, f'evector_J{j}_p{k}_i{i}.dat')
    #         for j in jrange for k in pars for i in self.nisotopes
    #     ]
    #     for efile in efiles:
    #         # nvib = np.arange(vrange[0], vrange[1] + 1, dtype=np.int64)
    #         evec = np.loadtxt(efile)
    #         arr = []
    #         for col in range(0, evec.shape[1]):
    #             for row in range(0, ninter):
    #                 index = 0
    #                 sum_arg = 0.0
    #                 for j in range(0, self.nch*self.ngrid):
    #                     if j % self.ngrid == 0:
    #                         index += 1
    #                     expr = igrid[row] - self.rgrid[j-index*self.ngrid]
    #                     if abs(expr) < 1e-15:
    #                         sinc = 1.0
    #                     else:
    #                         arg = (math.pi / istep) * expr
    #                         sinc = math.sin(arg) / arg
    #                     sum_arg += evec[j][col] * sinc
    #                 arr.append(sum_arg)
    #         wavefunc = np.array(arr).reshape(evec.shape[1], ninter).T
    #         nwavefunc = wavefunc / wavefunc.sum(axis=0, keepdims=1)
    #         self.save_wavefunctions(np.column_stack((igrid, nwavefunc)))

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

    def _save_output_data(self, stats=None):

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

        outd = np.c_[
            np.arange(1, self.out_data.shape[0]+1),
            self.out_data
        ]

        fmt = [
            '%7d', '%7d', '%7.1f', '%7.1f', '%7.1f', '%6d', '%6d',
            '%6d', '%15.6f', '%14.6f', '%12.6f', '%10.4f'
        ] + self.nch*['%8.3f'] + ['%5d']

        footer = ''
        if stats is not None:
            footer = \
                f'Chi Square = {stats[0]:.8f}\n' \
                f'RMS = {stats[1]:.8f} cm-1\n' \
                f'RMSD = {stats[2]:.8f}'

        with open(self.file_evals, 'w') as outs:
            np.savetxt(
                outs,
                outd,
                comments='#',
                header=header,
                footer=footer,
                fmt=fmt
            )

    def get_predicted_data(self):
        """Get all predicted eigenvalues as array

        Returns:
            numpy.ndarray: the computed eigenvalues
        """

        return self.evals_predicted

    def get_output_data(self):

        return self.out_data

    def get_Hamiltonian_matrix(self):

        # will return the matrix for the last comupted J, e/f level and isotope
        return self.hmatrix

    def save_predicted(self):

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

    def sort_predicted(self, cols=[]):

        out = []
        for col in cols:
            out.append(self.out_data[:, col])
        out = tuple(out)

        self.evals_predicted = self.evals_predicted[np.lexsort(out)]
        self.save_predicted()

    def sort_output(self, cols=[]):

        out = []
        for col in cols:
            out.append(self.out_data[:, col])
        out = tuple(out)

        self.out_data = self.out_data[np.lexsort(out)]
        self._save_output_data()

    def _group_increasing_sequences(self, iseq):

        return [
            list(arr) for arr in np.split(
                iseq, np.where(np.diff(iseq) != 1)[0] + 1
            ) if arr.size > 0
        ]

    def _find_divisors(self, numb):

        for i in [7, 6, 5, 4, 3, 2]:
            if numb % i == 0:
                return i
        return 1

    def save_data_info(self):

        s = ' '
        with open(self.info_outfile, 'w') as outf:
            outf.write(f'{30*"#"} Grid and Molecule Data {30*"#"}\n\n')
            outf.write(f'{s:>4}Number of Grid Points = {self.ngrid:<5d}\n\n')
            outf.write(
                f'{s:>4}Rmin = {self.rmin*C_bohr:>12.10f} '
                f'Angstrom = {self.rmin:>12.10f} Bohr\n'
            )
            outf.write(
                f'{s:>4}Rmax = {self.rmax*C_bohr:>12.10f} '
                f'Angstrom = {self.rmax:>12.10f} Bohr\n\n'
            )

            # rgrid_str = np.array2string(self.rgrid, precision=12)
            outf.write(f'{s:>4}Grid Points:\n\n')

            nd = self._find_divisors(self.rgrid.shape[0])
            np.savetxt(outf, self.rgrid.reshape(-1, nd), fmt='%20.14f')

            outf.write(f'\n{s:>4}Method of solution: {self.solver}\n\n')

            outf.write(
                f'{s:>4}Number of Defined Isotopes = {len(self.masses):>9d}\n'
            )
            outf.write(
                f'{s:>4}Number of Used Isotopes = {len(self.nisotopes):>12d}\n'
            )

            for n, niso in enumerate(self.nisotopes):
                mass = self.masses[n+1]
                outf.write(
                    f'{s:>6}Isotope {niso}: Mass = {mass:>15.9} au = '
                    f'{mass / C_massau:>15.9} amu\n'
                )

            outf.write(
                f'\n{s:>4}Number of Rotational Quanatum Numbers = '
                f'{self.jqnumbers.shape[0]}\n'
            )
            outf.write(f'\n{s:>4}Rotational Quanatum Numbers:\n')
            nd = self._find_divisors(self.jqnumbers.shape[0])
            np.savetxt(outf, self.jqnumbers.reshape(-1, nd), fmt='%7.1f')

            outf.write(f'\n{s:>4}Parity Levels = ')
            outf.write(', '.join(map(lambda x: 'e' if x else 'f', self.pars)))

            ofj = True if self.refj is not None else False
            outf.write(f'\n\n{s:>4}Shift Energies = {ofj}\n')
            if ofj:
                outf.write(f'{s:>4}Shift Energies by Level J = {self.refj}')

            outf.write(
                f'\n\n{s:>4}Hamiltonian Matrix Size = '
                f'{self.nch*self.ngrid}x{self.nch*self.ngrid}'
            )
            outf.write(
                f'\n{s:>4}Number of Calculated Eigenvalues = '
                f'{len(self.evals_predicted):>10d}'
            )
            outf.write(
                f'\n{s:>4}Number of Selected Eigenvalues = '
                f'{len(self.out_data):>12d}'
            )

            outf.write(f'\n\n{30*"#"} Channels Data {30*"#"}\n')
            outf.write(f'\n{s:>4}Number of Channels = {self.nch}\n')

            for ic, ch in enumerate(self.channels):
                outf.write(f'\n{s:>9} {ic+1}.\n')
                outf.write(f'{s:<13}Type: {ch.model}\n')
                outf.write(f'{s:<13}File: {ch.filep}\n')
                outf.write(f'{s:<13}Lambda: {ch.nlambda}\n')
                outf.write(f'{s:<13}Sigma: {ch.nsigma}\n')
                outf.write(f'{s:<13}Multiplicity: {ch.mult}\n')
                outf.write(f'{s:<13}Parameters:\n')
                if ch.model == 'pointwise':
                    for i in range(0, len(ch.R)):
                        outf.writelines(
                            f'{ch.R[i]:>29.14f}'
                            f'{ch.U[i]*C_hartree:>25.14f}\n'
                        )

            ugrid_cols = np.hstack((
                self.rgrid[:, np.newaxis] * Const.bohr,
                self.ugrid.reshape(self.nch, self.ngrid).T * C_hartree
            ))
            outf.write(f'\n{30*"#"} Channel Functions on Grid {30*"#"}\n\n')
            np.savetxt(outf, ugrid_cols, fmt=ugrid_cols.shape[1] * ['%21.12f'])

            outf.write(f'\n{30*"#"} Couplings Data {30*"#"}\n')
            outf.write(f'\n{s:>4}Number of Couplings = {self.ncp}\n')

            for ic, cp in enumerate(self.couplings):
                outf.write(f'\n{s:>9} {ic+1}.\n')
                outf.write(f'{s:<13}Type: {cp.model}\n')
                outf.write(f'{s:<13}Channels: {cp.interact}\n')
                outf.write(f'{s:<13}Coupling: {cp.coupling}\n')
                outf.write(f'{s:<13}Label: {cp.label}\n')
                outf.write(f'{s:<13}Multiplier: {cp.multiplier}\n')
                outf.write(f'{s:<13}Parameters:\n')
                if cp.model == 'pointwise':
                    for i in range(0, len(cp.xc)):
                        outf.writelines(
                            f'{cp.xc[i]:>29.14f}'
                            f'{cp.yc[i]:>25.14f}\n'
                        )

            fgrid_cols = np.hstack((
                self.rgrid[:, np.newaxis] * Const.bohr,
                self.fgrid.reshape(self.ncp, self.ngrid).T
            ))

            outf.write(f'\n{30*"#"} Coupling Functions on Grid {30*"#"}\n\n')
            np.savetxt(outf, fgrid_cols, fmt=fgrid_cols.shape[1] * ['%19.12f'])

            efmt = [
                '%7.1d', '%3.1d', '%7.1f', '%14.6f',
                '%4.1d', '%7.1f', '%4.1d', '%4.1d'
            ]

            outf.write(f'\n\n{20*"#"} Experimental data {20*"#"}\n')
            outf.write(
                f'\n{s:>4}File with Experimental Data = {self.exp_file}\n'
            )
            outf.write(f'\n{s:>4}Markers:')
            np.savetxt(outf, np.unique(self.exp_data[:, 6]), fmt='%3d')
            outf.write(
                f'\n{s:>4}Number of used experimental data = '
                f'{self.exp_data.shape[0]}\n\n'
            )

            np.savetxt(outf, self.exp_data, comments='', fmt=efmt)

            avg_unc = np.sum(self.exp_data[:, 5]) / self.exp_data.shape[0]
            outf.write(
                f'\n{s:>4}Average experimental uncertainty = '
                f'{avg_unc:7.4f}'
            )

            outf.write(f'\n\n{20*"#"} Eigenvalues {20*"#"}\n\n')
            outf.write(
                f'{s:>4}Number of computed eigenvalues = '
                f'{self.out_data.shape[0]}\n\n'
            )

            outd = np.c_[
                np.arange(1, self.out_data.shape[0]+1),
                self.out_data
            ]

            fmt = [
                '%7d', '%7d', '%7.1f', '%7.1f', '%7.1f', '%6d', '%6d',
                '%6d', '%15.6f', '%14.6f', '%12.6f', '%10.4f'
            ] + self.nch*['%8.3f'] + ['%5d']

            np.savetxt(outf, outd, fmt=fmt)

            chi2, rms, rmsd = self.calculate_stats(
                self.out_data[:, 8], self.out_data[:, 7],
                self.out_data[:, 10], False
            )

            outf.write(
                f'\n{s:>4}Chi Square = {chi2:.8f}\n'
                f'{s:>4}RMS = {rms:.8f} cm-1\n'
                f'{s:>4}RMSD = {rmsd:.8f}'
            )

            outf.write(f'\n\n{s:>4}Levels with largest deviations:\n\n')

            outf.write(f'\n\n{20*"#"} Predicted Eigenvalues {20*"#"}\n\n')
