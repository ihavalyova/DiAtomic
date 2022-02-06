import io
import numpy as np
import wavenumbers
import wavenumbers_intens

from scipy.interpolate import CubicSpline
from utils import C_bohr, C_boltzmannk
try:
    import py3nj
except ModuleNotFoundError:
    pass
    # print("'py3nj' module is not installed!\n")


class Spectrum:
    """[summary]
    """

    def __init__(self, H, dmfs=None, spec_type='absorption'):

        self.H = H
        self.rgrid = self.H.rgrid
        self.rmax = H.rmax
        self.rmin = H.rmin
        self.ngrid = H.ngrid
        self.rstep = H.rstep

        self.dmfs_init = None
        if dmfs is not None:
            self.dmfs_init = {}
            for (n, k) in dmfs:
                self.dmfs_init[(n, k)] = np.loadtxt(io.StringIO(dmfs[(n, k)]).getvalue())

        self.wavens = None
        self.hlf = None
        self.spec_type = spec_type

        # default file names
        self.fname_wavens = 'out_wavenumbers.dat'
        self.fname_hlf = 'hlf.dat'
        self.fname_lifetimes = 'lifetimes.dat'
        self.fname_rel_int = 'relative_intensities.dat'
        self.fname_acoefs = 'Acoefs.dat'
        self.fname_compare = 'compared_wavenumbers.dat'

    def set_constraints(self, uenr_range=None, lenr_range=None, lsymm=None,
                        usymm=None, ujrange=None, ljrange=None, wrange=None):

        # set default constraints
        self.freq_range = wrange or (0, 1e5)
        usymm = usymm or (0, 1)
        lsymm = lsymm or (0, 1)

        evalues = self.H.calc_data

        eind, jind, pind = 1, 2, 3

        # constraints by energy
        self.lower_levels = evalues[(evalues[:, eind] >= lenr_range[0]) &
                                    (evalues[:, eind] <= lenr_range[1])]

        self.upper_levels = evalues[(evalues[:, eind] >= uenr_range[0]) &
                                    (evalues[:, eind] <= uenr_range[1])]

        # constraints by J
        ujs, uje = ujrange[0], ujrange[1]

        # upper_jrots = np.arange(ujs, uje+1)
        self.upper_levels = self.upper_levels[(self.upper_levels[:, jind] >= ujs) &
                                              (self.upper_levels[:, jind] <= uje)]

        ljs, lje = ljrange[0], ljrange[1]
        # lower_jrots = np.arange(ljs, lje+1)
        self.lower_levels = self.lower_levels[
            (self.lower_levels[:, jind] >= ljs) &
            (self.lower_levels[:, jind] <= lje)
        ]

        # constarints by symmetry
        self.lower_levels = self.lower_levels[
            (np.in1d(self.lower_levels[:, pind], lsymm[0])) |
            (np.in1d(self.lower_levels[:, pind], lsymm[1]))
        ]
        self.upper_levels = self.upper_levels[
            (np.in1d(self.upper_levels[:, pind], usymm[0])) |
            (np.in1d(self.upper_levels[:, pind], usymm[1]))
        ]

    def calculate_wavenumbers(self, ulevels, llevels, wrange=None, apply_rules=False):

        if self.H.wavens_data is None:
            return self._calculate_wavenumbers_without_observations(
                ulevels, llevels, wrange, apply_rules)
        else:
            return self._calculate_wavenumber_with_observations(
                ulevels, llevels, wrange, apply_rules)

    def _calculate_wavenumbers_without_observations(self, ulevels, llevels, wrange, apply_rules):

        # id v J E p s iso lambda omega
        ulevels = ulevels[:, [0, -3, 2, 1, 3, -4, 4, -2, -1]]
        llevels = llevels[:, [0, -3, 2, 1, 3, -4, 4, -2, -1]]

        # convert to row-major/c-style/c-contiguous format
        ulevels = np.ascontiguousarray(ulevels)
        llevels = np.ascontiguousarray(llevels)

        # uid uv uJ up us ucE ul uo lid lv lJ lp ls lcE ll lo cw branch
        out_wavens = wavenumbers.calculate_wavenumbers_list(ulevels, llevels)

        # apply strict selection rules by J and symmetry through the branch label
        out_wavens = out_wavens[out_wavens[:, -1] != -1]

        # filter by wavenumber value
        if wrange is not None:
            out_wavens = out_wavens[out_wavens[:, -2] >= wrange[0] &
                                    out_wavens[:, -2] <= wrange[1]]

        # apply non-strict selection rules by computing the Honl-London factors
        if apply_rules:
            self.hlf = self._apply_rules(out_wavens)
            out_wavens = out_wavens[np.where(self.hlf != 0.0)[0], :]
            self.hlf = self.hlf[np.where(self.hlf != 0.0)[0]]

        return out_wavens

    def _calculate_wavenumber_with_observations(self, ulevels, llevels, wrange=None,
                                                apply_rules=False):

        # id v J E p s iso lambda omega
        ulevels = ulevels[:, [0, -3, 2, 1, 3, -4, 4, -2, -1]]
        llevels = llevels[:, [0, -3, 2, 1, 3, -4, 4, -2, -1]]

        # find the matching experimental and calculated term values
        euterms = np.unique(self.H.wavens_data[:, :4], axis=0)
        elterms = np.unique(self.H.wavens_data[:, 4:8], axis=0)
        v, j, pi, si = 1, 2, 4, 5
        uinds, linds = [], []
        for i in euterms:
            inds = np.where(np.all(ulevels[:, [v, j, pi, si]] == i, axis=1))[0]
            if len(inds) > 0:
                uinds.append(inds[0])
        for i in elterms:
            inds = np.where(np.all(llevels[:, [v, j, pi, si]] == i, axis=1))[0]
            if len(inds) > 0:
                linds.append(inds[0])
        llevels = llevels[linds, :]
        ulevels = ulevels[uinds, :]

        # convert to row-major/c-style/c-contiguous format
        ulevels = np.ascontiguousarray(ulevels)
        llevels = np.ascontiguousarray(llevels)

        # uid uv uJ up us ucE ul uo lid lv lJ lp ls lcE ll lo cw ew dw uncw eint unci branch
        out_wavens = wavenumbers_intens.calculate_wavenumbers_list(
            ulevels, llevels, self.H.wavens_data)

        # apply strict selection rules by J and symmetry through the branch label
        out_wavens = out_wavens[out_wavens[:, -1] != -1]

        # chi2, rms, rmsd = self.H.calculate_stats(out_wavens[:, 17], out_wavens[:, 16],
        #                                          out_wavens[:, 19], is_weighted=False)

        # apply non-strict selection rules by computing the Honl-London factors
        if apply_rules:
            self.hlf = self._apply_rules(out_wavens)
            out_wavens = out_wavens[np.where(self.hlf != 0.0)[0], :]
            self.hlf = self.hlf[np.where(self.hlf != 0.0)[0]]

        return out_wavens

    def wrapper_calculate_wavenumbers(self, ypar):

        self.H.interpolate_functions(ypar)
        self.H.solve(energy_subset_index=self.H.energy_subset_index,
                     energy_subset_value=self.H.energy_subset_value)
        ulevels, llevels = self.H.extract_terms_in_range(
            uenergy=self.H.uenergy, lenergy=self.H.lenergy, usymm=self.H.usymm,
            lsymm=self.H.lsymm, uj=self.H.uj, lj=self.H.lj,
            ustate=self.H.ustate, lstate=self.H.lstate)

        return self._calculate_wavenumber_with_observations(ulevels, llevels)

    def _apply_rules(self, out_wavens):

        usymm, lsymm = out_wavens[:, 3], out_wavens[:, 11]
        ujq, ljq = out_wavens[:, 2], out_wavens[:, 10]
        ulambda, llambda = out_wavens[:, 6], out_wavens[:, 14]
        uomega, lomega = out_wavens[:, 7], out_wavens[:, 15]

        return self._compute_honl_london_factors(
            usymm, lsymm, ujq, ljq, uomega, lomega, ulambda, llambda)

    def calculate_Einstein_coefficients(self, out_wavens, dmfs=None, ninter=1000):
        # if self.hlf is None:
        #     self._apply_rules(out_wavens)

        if dmfs is None:
            dmfs = self.dmfs_init

        # uid uv uJ up us ucE ul uo lid lv lJ lp ls lcE ll lo cw ew dw uncw eint unci branch
        self.edipole_element = self._compute_line_strength(out_wavens, dmfs, ninter=ninter)

        # np.savetxt('dipole_element.dat', self.edipole_element, fmt='%14.8e')
        line_strength = self.edipole_element[:, -1] * self.hlf
        jinitial = self.edipole_element[:, 2]
        waven = self.edipole_element[:, 16]
        # statistical weight of the initial level
        self.gji = (2 * jinitial + 1)
        self.acoef = self._calculate_Einstein_coeffcients(waven, line_strength, self.gji)
        acoef_result = np.c_[self.edipole_element[:, :-1], self.acoef]
        self.nonzero_ainds = np.where(self.acoef != 0.0)[0]
        self.acoef_final = acoef_result[self.nonzero_ainds, :]

        return self.acoef_final

    def wrapper_calculate_Einstein_coefficients(self, ypar):

        # TODO: call wrapper_calculate_wavenumbers()
        self.H.interpolate_functions(ypar)
        self.H.solve(energy_subset_index=self.H.energy_subset_index,
                     energy_subset_value=self.H.energy_subset_value)
        ulevels, llevels = self.H.extract_terms_in_range(
            uenergy=self.H.uenergy, lenergy=self.H.lenergy, usymm=self.H.usymm,
            lsymm=self.H.lsymm, uj=self.H.uj, lj=self.H.lj,
            ustate=self.H.ustate, lstate=self.H.lstate)
        out_wavens = self.calculate_wavenumbers(ulevels, llevels)

        dmfs = None
        if self.dmfs_init is not None:
            dmfs = {}
            ypar_shape = self.H.ypar_init.shape[0]
            for (n, k) in self.dmfs_init:
                dmf_params = ypar[ypar_shape:self.dmfs_init[(n, k)].shape[0]+ypar_shape]
                dmfs[(n, k)] = np.column_stack((self.dmfs_init[(n, k)][:, 0], dmf_params,
                                                self.dmfs_init[(n, k)][:, 2]))

        return self.calculate_Einstein_coefficients(out_wavens, dmfs=dmfs)

    def _calculate_Einstein_coeffcients(self, wavenumber, line_strength, gji):

        # e0 = constants.value('vacuum electric permittivity')
        # planck = constants.value('Planck constant')
        # consts = (16 * np.pi**3) / 3 * e0 * planck
        # acoef = (consts * wavenumber**3 * line_strength) / gi
        acoef = (3.1361891e-7 * line_strength * wavenumber**3) / gji

        return acoef

    def _compute_line_strength(self, out_wavens, dmfs, ninter=1000):

        ivec_inds = out_wavens[:, 0].astype(np.int)
        fvec_inds = out_wavens[:, 8].astype(np.int)

        rgrid, rstep = self.rgrid * C_bohr, self.rstep * C_bohr
        igrid, istep = np.linspace(self.rmin, self.rmax, num=ninter, endpoint=True, retstep=True)
        igrid, istep = igrid * C_bohr, istep * C_bohr

        sinc_matrix = self._calculate_sinc_matrix(rgrid, igrid, rstep)
        dipole_matrix = self._interpolate_dipole_moment(igrid, dmfs)

        # np.savetxt('sinc_matrix.dat', sinc_matrix, fmt='%14.11e')
        # np.savetxt('evecs_all.dat', self.H.evecs_matrix, fmt='%12.11e')
        # np.savetxt('dipole_matrix.dat', dipole_matrix[0, 1, :], fmt='%14.8e')
        result = np.zeros((ivec_inds.shape[0], out_wavens.shape[1]+1))

        ii = 0
        for (i, j) in zip(ivec_inds, fvec_inds):
            dme = 0.0
            for n in range(self.H.nch):
                for m in range(self.H.nch):
                    ncoefs = self.H.evecs_matrix[n*self.H.ngrid:(n+1)*self.H.ngrid, i-1]
                    kcoefs = self.H.evecs_matrix[m*self.H.ngrid:(m+1)*self.H.ngrid, j-1]
                    # for l in range(ngrid):
                    #     for p in range(ngrid):
                    #         sincl = sinc_matrix[l, :]
                    #         sincp = sinc_matrix[p, :]
                    #         sumq = np.sum(sincl*dipole_matrix[n, m, :]*sincp)
                    #         dme += kcoefs[l] * ncoefs[p] * sumq * istep
                    # res = dme / rstep
                    # res = dme / rstep # this should be outside n, m loops
                    sumq = np.dot(sinc_matrix, (dipole_matrix[n, m, :]*sinc_matrix).T)
                    dme += np.dot(kcoefs, ncoefs * sumq) * istep
            res = np.sum(dme) / rstep
            result[ii, :] = np.hstack((out_wavens[ii, :], res**2))
            ii += 1

        return result

    def _interpolate_dipole_moment(self, igrid, dmfs):

        #dipole_matrix = np.zeros((self.H.nch, self.H.nch, igrid.shape[0]))
        dipole_matrix = np.ones((self.H.nch, self.H.nch, igrid.shape[0]))

        for n in range(1, self.H.nch+1):
            for k in range(1, self.H.nch+1):
                is_interp = False
                try:
                    dmx, dmy = dmfs[(n,k)][:, 0], dmfs[(n,k)][:, 1]
                    cs = CubicSpline(dmx, dmy, bc_type='natural')
                    dmi = cs(igrid)
                    is_interp = True
                except (KeyError, TypeError):
                    # dmi = np.ones_like(igrid)
                    continue

                dipole_matrix[n-1, k-1, :] = dmi
                if n != k and is_interp:
                    dipole_matrix[k-1, n-1, :] = dmi

        return dipole_matrix

    def _calculate_sinc_matrix(self, rgrid, igrid, rstep):

        sinc_matrix = np.zeros((self.H.ngrid, igrid.shape[0]))

        for i in range(self.H.ngrid):
            argi = (igrid - rgrid[i]) / rstep
            sinc = np.sinc(argi)
            sinc_matrix[i, :] = sinc

        return sinc_matrix

    def calculate_relative_intensities(self, wavens, T=296):

        if self.spec_type.lower().startswith('a'):

            # TODO: check if Einstein coefs have already been calculated
            kt = C_boltzmannk * T

            # TODO check column numbers
            intensity = self.acoef * self.gji * wavens[:, -1] * np.exp(-wavens[:, 10]/kt)
            # (1 - np.exp(-freq_calc[:, -1]/kt))
            # np.square(np.square(freq_calc[:, -1]))
            intensity_result = np.c_[self.edipole_element[:, :-1], intensity]
            intensity_final = intensity_result[self.nonzero_ainds, :]
        else:
            pass

        return intensity_final

    def calculate_lifetimes(self, save=False, filename=None):

        # t = 1 / sum_{j} A_ij

        # get the indices of the unique upper levels
        _, unq_uinds, unq_uinv = np.unique(
            self.acoef_final[:, :6], return_index=True,
            return_inverse=True, axis=0)

        # sum the Einstein coefficients for each group of unique upper levels
        # sum_acoefs = np.zeros(unq_uinv.shape[0])
        sum_acoefs = np.zeros(np.unique(unq_uinv).shape)

        for i, ii in enumerate(unq_uinv):
            sum_acoefs[ii] += self.acoef_final[i, -1]

        # the lifetimes are the inverse of the acumulated sums
        lifetimes = 1.0 / sum_acoefs

        # conncatenate the upper levels and the calculated lifetimes
        lifetimes_final = np.column_stack((
            self.acoef_final[unq_uinds, :6],
            lifetimes[:unq_uinds.shape[0]]))

        if save:
            filename = filename or self.fname_lifetimes
            self._save_lifetimes(lifetimes_final, filename)

        return lifetimes_final

    def calculate_partition_functions(self, Ia=0, Ib=0, T=296,
                                      save=False, filename=None):
        # for calculating the partition functions the total degeneracy factor
        # for the molecule is needed which includes also the nuclear degeneracy

        # Q(T) = g_ns sum_{i} (2J_{i}+1) * exp(-E_{i}/kT)
        # for heteronuclear diatomics: g_ns = (2I_{a} + 1)(2I_{b} + 1)
        # for homonuclear diatomics: depends on the nuclear symmetry

        gns = 1
        # Ia = Ib => the molecule is homonuclear
        if Ia == Ib:
            # if Ia.is_integer():
            raise NotImplementedError(
                'Homonuclear diatomics are currently not supported')
        else:
            gns = (2*Ia+1)*(2*Ib+1)

        kt = C_boltzmannk * T
        # qpart = gns *

    def calculate_absolute_intensity(self):

        # the formula for the absolute line intensity is diffrent for
        # absorption and emission

        if self.spectrum == 'absorption' or self.spectrum == 'a':
            self._calculate_absorption_line_intensity()
        else:
            self._calculate_emission_line_intensity()

    def _calculate_absorption_line_intensity(self):
        pass

    def _calculate_emission_line_intensity(self):
        pass

    def calculate_oscilator_strength(self):
        raise NotImplementedError()

    def _compute_honl_london_factors(self, usymm, lsymm, uj, lj, uomega,
                                     lomega, ulambda, llambda):

        # n' E' J' p' iso' st' v' l' o' n E J p iso st v l o freq

        ueps = self._map_parity_to_epsilon(usymm)
        leps = self._map_parity_to_epsilon(lsymm)
        eps_expr = 1.0 + (ueps * leps * np.power(-1, 1.0 + uj - lj))

        udelta = self._map_quantum_number_to_delta(ulambda)
        udelta *= self._map_quantum_number_to_delta(uomega)
        ldelta = self._map_quantum_number_to_delta(llambda)
        ldelta *= self._map_quantum_number_to_delta(lomega)
        delta_expr = 1.0 + udelta + ldelta - 2.0 * udelta * ldelta

        j_expr = ((2.0 * uj) + 1.0) * ((2.0 * lj + 1.0))

        two_l1 = np.int64(2*uj)
        two_l2 = 2*np.ones(uj.shape[0], dtype=np.int64)
        two_l3 = np.int64(2*lj)
        two_m1 = np.int64(-2*uomega)
        two_m2 = np.int64(2*(ulambda-llambda))
        two_m3 = np.int64(2*lomega)

        # allows the ambiguous sign in the 3j symbol when one of the Lambda
        # quantum numbers is zero (Ref: J.K.G. Watson, JMS 252 (2008))

        if (ulambda.any() == 0 and llambda.any() != 0) or \
           (ulambda.any() != 0 and llambda.any() == 0):
            two_m2 = np.int64(2*(ulambda+llambda))
            two_m3 = np.int64(-2*lomega)

        # qunatum numbers that do not satisfy the following
        # conditions should be set to zero

        if (two_l1 < np.abs(two_m1)).any():
            valid = 1*~(two_l1 < np.abs(two_m1))
            two_m1 *= valid

        if (two_l2 < np.abs(two_m2)).any():
            valid = 1*~(two_l2 < np.abs(two_m2))
            two_m2 *= valid

        if (two_l3 < np.abs(two_m3)).any():
            valid = 1*~(two_l3 < np.abs(two_m3))
            two_m3 *= valid

        wigner_3j = py3nj.wigner3j(
            two_l1, two_l2, two_l3, two_m1, two_m2, two_m3
        )
        wigner_3j_square = wigner_3j ** 2

        hlf = 0.5 * eps_expr * delta_expr * j_expr * wigner_3j_square

        return hlf

    def _map_parity_to_epsilon(self, pars):

        return np.array(list(map(lambda x: -1 if x == 0 else x, pars)))

    def _map_quantum_number_to_delta(self, qlist):

        return np.array(list(map(lambda x: 1 if x == 0 else 0, qlist)))

    def save_wavenumbers(self, out_wavens, filename=None):

        if self.H.wavens_data is None:
            self._save_wavenumbers_without_observations(out_wavens, filename)
        else:
            self._save_wavenumbers_with_obseravtions(out_wavens, filename)

    def _save_wavenumbers_with_obseravtions(self, out_wavens, filename):

        filename = filename or self.fname_wavens

        cols = list(range(1,8)) + list(range(9, out_wavens.shape[1]))
        wavens = out_wavens[:, cols]

        labels = ("v'", "J'", "symm'", "state'", "E'", "Lambda'", "Omega'",
                  'v', 'J', 'symm', 'state', 'E', 'Lambda', 'Omega', 'cfreq',
                  'efreq', 'diff_freq', 'unc_freq', 'eintens', 'unc_int', 'branch')

        header = (f'{labels[0]:^7}{labels[1]:<5}{labels[2]:<5}{labels[3]:<10}'
                  f'{labels[4]:<6}{labels[5]:<8}{labels[6]:<7}{labels[7]:<6}'
                  f'{labels[8]:<4}{labels[9]:<5}{labels[10]:<10}{labels[11]:<6}'
                  f'{labels[12]:<7}{labels[13]:<8}{labels[14]:<10}{labels[15]:<9}'
                  f'{labels[16]:<10}{labels[17]:<14}{labels[18]:<11}{labels[19]:<10}'
                  f'{labels[20]}')

        fmt = ('%5.1d', '%6.1f',  '%4.1d', '%4.1d', '%12.5f', '%4.1d',
               '%5.1f', '%4.1d', '%5.1f', '%4.1d', '%4.1d', '%12.5f',
               '%4.1d', '%4.1f', '%15.8f', '%15.8f', '%11.3e', '%7.3f',
               '%14.5e', '%8.2f', '%6.1d')

        np.savetxt(filename, wavens, header=header, fmt=fmt)

    def _save_wavenumbers_without_observations(self, out_wavens, filename):

        filename = filename or self.fname_wavens
        cols = list(range(1,8)) + list(range(9, out_wavens.shape[1]))
        wavens = out_wavens[:, cols]

        labels = ("v'", "J'", "symm'", "state'", "E'", "Lambda'", "Omega'",
                  "v", "J", "symm", "state", "E", "Lambda", "Omega", "cfreq", "branch")

        header = (f'{labels[0]:^7}{labels[1]:<5}{labels[2]:<5}{labels[3]:<10}'
                  f'{labels[4]:<6}{labels[5]:<8}{labels[6]:<7}{labels[7]:<6}'
                  f'{labels[8]:<4}{labels[9]:<5}{labels[10]:<10}{labels[11]:<6}'
                  f'{labels[12]:<7}{labels[13]:<8}{labels[14]:<8}{labels[15]}')

        fmt = ('%5.1d', '%6.1f',  '%4.1d', '%4.1d', '%12.5f', '%4.1d',
               '%5.1f', '%4.1d', '%5.1f', '%4.1d', '%4.1d', '%12.5f',
               '%4.1d', '%4.1f', '%12.5f', '%4.1d')

        np.savetxt(filename, wavens, header=header, fmt=fmt)

    def save_Honl_London_factor(self, out_wavens, hlf, filename=None):

        calc_hlf = np.c_[out_wavens, hlf]
        # calc_hlf = calc_hlf[np.where(hlf != 0.0)[0], :]
        # cols = [6, 2, 1, 3, 4, 5, 15, 11, 10, 12, 13, 14, -2, -1]
        filename = filename or self.fname_hlf

        labels = ("v'", "J'", "E'", "symm'", "iso'", "state'", 'v',
                  'J', 'E', 'symm', 'iso', 'state', 'freq', 'hlf')

        header = (f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<6}'
                  f'{labels[4]:<6}{labels[5]:<9}{labels[6]:<6}{labels[7]:<13}'
                  f'{labels[8]:<9}{labels[9]:<6}{labels[10]:<6}'
                  f'{labels[11]:<10}{labels[12]:<14}{labels[13]}')

        fmt = ('%6.1d', '%7.1f', '%14.5f', '%5.1d', '%5.1d',
               '%5.1d', '%7.1d', '%7.1f', '%14.5f', '%5.1d',
               '%5.1d', '%5.1d', '%15.5f', '%12.5f')

        np.savetxt(filename, calc_hlf, header=header, fmt=fmt)

    def save_rel_intensity(self, result, filename=None):

        filename = filename or self.fname_rel_int

        labels = ("v'", "J'", "E'", "symm'", "iso'", "state'", 'v',
                  'J', 'E', 'symm', 'iso', 'state', 'freq', 'intensity')

        header = (f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<6}'
                  f'{labels[4]:<6}{labels[5]:<9}{labels[6]:<6}{labels[7]:<13}'
                  f'{labels[8]:<9}{labels[9]:<6}{labels[10]:<6}'
                  f'{labels[11]:<10}{labels[12]:<15}{labels[13]}')

        fmt = ('%5.1d', '%7.1f', '%14.5f', '%5.1d', '%5.1d',
               '%5.1d', '%7.1d', '%7.1f', '%14.5f', '%5.1d',
               '%5.1d', '%5.1d', '%15.5f', '%16.5e')

        np.savetxt(filename, result, comments='#', header=header, fmt=fmt)

    def save_einstein_coeffcients(self, acoefs, filename=None):

        if self.H.wavens_data is None:
            self._save_einstein_coeffcients_without_observations(acoefs, filename)
        else:
            self._save_einstein_coeffcients_with_observations(acoefs, filename)

    def _save_einstein_coeffcients_without_observations(self, acoefs, filename):

        filename = filename or self.fname_acoefs

        cols = list(range(1,8)) + list(range(9, acoefs.shape[1]))
        acoefs_out = acoefs[:, cols]

        labels = ("v'", "J'", "symm'", "state'", "E'", "Lambda'", "Omega'", "v",
                  "J", "symm", "state", "E", "Lambda", "Omega", "cfreq", "branch", "A")

        header = (f'{labels[0]:^7}{labels[1]:<6}{labels[2]:<6}{labels[3]:<10}'
                  f'{labels[4]:<7}{labels[5]:<8}{labels[6]:<7}{labels[7]:<6}'
                  f'{labels[8]:<6}{labels[9]:<6}{labels[10]:<10}{labels[11]:<8}'
                  f'{labels[12]:<7}{labels[13]:<8}{labels[14]:<14}{labels[15]}')

        fmt = ('%5.1d', '%6.1f',  '%5.1d', '%5.1d', '%12.5f', '%5.1d',
               '%5.1f', '%5.1d', '%5.1f', '%5.1d', '%5.1d', '%12.5f',
               '%5.1d', '%5.1f', '%12.5f', '%8.1d', '%12.5e')

        np.savetxt(filename, acoefs_out, comments='#', header=header, fmt=fmt)

    def _save_einstein_coeffcients_with_observations(self, acoefs, filename):

        filename = filename or self.fname_acoefs

        # the diff column for A is added here
        adiff = acoefs[:, -4] - acoefs[:, -1]
        acoefs_out = np.column_stack((acoefs, adiff))

        cols = list(range(1,8)) + list(range(9, acoefs_out.shape[1]))
        acoefs_out = acoefs_out[:, cols]

        labels = ("v'", "J'", "symm'", "state'", "E'", "Lambda'", "Omega'",
                  "v", "J", "symm", "state", "E", "Lambda", "Omega", "cfreq", "efreq",
                  "diff_freq", "unc_freq", "eintens", "unc_int", "branch", "A", "diff_A")

        header = (f'{labels[0]:^7}{labels[1]:<5}{labels[2]:<5}{labels[3]:<10}'
                  f'{labels[4]:<6}{labels[5]:<8}{labels[6]:<7}{labels[7]:<6}'
                  f'{labels[8]:<4}{labels[9]:<5}{labels[10]:<10}{labels[11]:<6}'
                  f'{labels[12]:<7}{labels[13]:<8}{labels[14]:<12}{labels[15]:<10}'
                  f'{labels[16]:<10}{labels[17]:<10}{labels[18]:<11}{labels[19]:<9}'
                  f'{labels[20]:<11}{labels[21]:<10}{labels[22]}')

        fmt = ('%5.1d', '%6.1f',  '%4.1d', '%4.1d', '%12.5f', '%4.1d',
               '%5.1f', '%4.1d', '%5.1f', '%4.1d', '%4.1d', '%12.5f',
               '%4.1d', '%4.1f', '%15.8f', '%15.8f', '%11.3e', '%7.3f',
               '%16.10e', '%8.2f', '%5.1d', '%16.10e', '%12.3e')

        np.savetxt(filename, acoefs_out, comments='#', header=header, fmt=fmt)

    def save_lifetimes(self, lifetimes, filename=None):

        labels = ("v'", "J'", "E'", "symm'", "iso'", "state'", 'life_time')

        header = (f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<8}'
                  f'{labels[4]:<7}{labels[5]:<10}{labels[6]}')

        fmt = ('%5.1d', '%7.1f', '%14.5f', '%6.1d', '%6.1d', '%6.1d', '%17.5e')

        np.savetxt(filename, lifetimes, comments='#', header=header, fmt=fmt)

    def _get_default_jrange(self, uterms, uj, lterms, lj):

        js = min(np.amin(uterms[:, uj]), np.amin(lterms[:, lj]))
        je = max(np.amax(uterms[:, uj]), np.amax(lterms[:, lj]))
        return js, je

    def _map_number_to_branch(self):

        return {
            1: 'Pee',
            2: 'Pff',
            3: 'Qef',
            4: 'Qfe',
            5: 'Ree',
            6: 'Rff'
        }

    def _map_branch_to_number(self):

        return {
            'Pee': 1,
            'Pff': 2,
            'Qef': 3,
            'Qfe': 4,
            'Ree': 5,
            'Rff': 6
        }
