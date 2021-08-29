import io
import numpy as np
from scipy.interpolate import CubicSpline
from utils import C_bohr
from scipy import constants

try:
    import py3nj
except ModuleNotFoundError:
    pass
    # print("'py3nj' module is not installed!\n")


class Spectrum:

    def __init__(self, mlevels, grid, spectrum='absorption'):

        # TODO: allow for setting an input file with computed energies and grid
        self.mlevels = mlevels
        self.grid = grid

        self.dfreq_range = (0, 1e5)
        self.dusymm = (0, 1)
        self.dlsymm = (0, 1)
        self.is_hlf_computed = False
        self.freq_calc = None
        self.spectrum = spectrum
        # TODO: set default constarints for J, par, El, Eu, freq

        # default file names
        self.fname_wavenumbers = 'wavenumbers.dat'
        self.fname_hlf = 'hlf.dat'
        self.fname_lifetimes = 'lifetimes.dat'
        self.fname_rel_int = 'relative_intensities.dat'
        self.fname_acoefs = 'Acoefs.dat'
        self.fname_compare = 'compared_wavenumbers.dat'

        self.c_boltzmannk = constants.value(
            'Boltzmann constant in inverse meter per kelvin'
        ) * 0.01

    def set_constraints(self, Elower=None, Eupper=None, lsymm=None, usymm=None,
                        ujrange=None, ljrange=None, freq_range=None):

        # set default constraints
        self.freq_range = freq_range or (0, 1e5)
        usymm = usymm or self.dusymm
        lsymm = lsymm or self.dlsymm

        # TODO: set default constarints for J, par, El, Eu, freq

        evalues = self.mlevels.get_predicted_data()
        eind, jind, pind = 1, 2, 3

        # constraints by energy
        self.lower_levels = evalues[
            (evalues[:, eind] >= Elower[0]) & (evalues[:, eind] <= Elower[1])
        ]
        self.upper_levels = evalues[
            (evalues[:, eind] >= Eupper[0]) & (evalues[:, eind] <= Eupper[1])
        ]

        # constraints by J
        ujs, uje = ujrange[0], ujrange[1]

        # upper_jrots = np.arange(ujs, uje+1)
        self.upper_levels = self.upper_levels[
            (self.upper_levels[:, jind] >= ujs) &
            (self.upper_levels[:, jind] <= uje)
        ]

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

    def calculate_wavenumbers(self, apply_rules=False, array_size=10000,
                              save=False, filename=None):

        jind, pind = 2, 3
        # TODO: remove self parameters from the function call

        # create a list of wavenumbers based on the strict selection rules
        # for J and symmetry
        self.freq_calc = self._create_wavenumbers_list(
            self.upper_levels, self.lower_levels, jind,
            pind, self.freq_range, array_size
        )

        # apply additional selection rules by computing Honl-London factors
        if apply_rules:
            self._apply_rules_by_HLF()

        cols = [6, 2, 1, 3, 4, 5, 15, 11, 10, 12, 13, 14, -1]
        if save:
            filename = filename or self.fname_wavenumbers
            self._save_wavenumbers(self.freq_calc[:, cols], filename)

        return self.freq_calc[:, cols]

    def _apply_rules_by_HLF(self):

        self.hlf = np.ones(self.freq_calc.shape[0])

        # extract the quantum numbers from the computed wavenumbers list
        usymm, lsymm = self.freq_calc[:, 3], self.freq_calc[:, 12]
        uj, lj = self.freq_calc[:, 2], self.freq_calc[:, 11]
        ulambda, llambda = self.freq_calc[:, 7], self.freq_calc[:, -3]
        uomega, lomega = self.freq_calc[:, 8], self.freq_calc[:, -2]

        self.hlf = self._compute_Honl_London_factors(
            usymm, lsymm, uj, lj, uomega, lomega, ulambda, llambda
        )

        self.is_hlf_computed = True

    def save_HLF(self, filename=None):

        calc_hlf = np.c_[self.freq_calc, self.hlf]
        # calc_hlf = calc_hlf[np.where(hlf != 0.0)[0], :]
        cols = [6, 2, 1, 3, 4, 5, 15, 11, 10, 12, 13, 14, -2, -1]
        filename = filename or self.fname_hlf
        self._save_Honl_London_factor(calc_hlf[:, cols], filename)

    def calculate_Einstein_coefficients(self, ninter=1000, dmf=None,
                                        save=False, filename=None):

        # TODO: check if HLF has already been computed or not
        # TODO: check if freq_calc has already been computed or not

        initial_levels = self.freq_calc[:, :9]
        final_levels = self.freq_calc[:, 9:-1]

        self.edipole_element = self._compute_electric_dipole_elements(
            initial_levels, final_levels, dmf, ninter=ninter
        )

        line_strength = self.edipole_element[:, -1] * self.hlf

        jinitial = self.edipole_element[:, 1]
        wavenumber = self.edipole_element[:, -2]

        # statistical weight of the initial level
        self.gji = (2 * jinitial + 1)
        self.acoef = self._calculate_Einstein_coeffcients(
            wavenumber, line_strength, self.gji
        )
        acoef_result = np.c_[self.edipole_element[:, :-1], self.acoef]
        self.nonzero_ainds = np.where(self.acoef != 0.0)[0]
        self.acoef_final = acoef_result[self.nonzero_ainds, :]

        if save:
            filename = filename or self.fname_acoefs
            self._save_Einstein_coeffcients(self.acoef_final, filename)

        return self.acoef_final

    def calculate_relative_intensities(self, T=296, save=False, filename=None):

        if self.spectrum == 'absorption' or self.spectrum == 'a':

            # TODO: check if Einstein coefs have already been calculated

            kt = self.c_boltzmannk * T
            print(self.freq_calc[:, 10])
            print(self.freq_calc[:, -1])
            # TODO check this
            intensity = self.acoef * self.gji * self.freq_calc[:, -1] * np.exp(-self.freq_calc[:, 10]/kt)
            # (1 - np.exp(-freq_calc[:, -1]/kt))
            # np.square(np.square(freq_calc[:, -1]))
            intensity_result = np.c_[self.edipole_element[:, :-1], intensity]
            intensity_final = intensity_result[self.nonzero_ainds, :]

            # cols = [6, 2, 1, 3, 4, 5, 15, 11, 10, 12, 13, 14, -2, -1]
            if save:
                filename = filename or self.fname_rel_int
                self._save_rel_intensity(intensity_final, filename)
        else:
            pass

        return intensity_final

    def calculate_lifetimes(self, save=False, filename=None):

        # t = 1 / sum_{j} A_ij

        # get the indices of the unique upper levels
        _, unq_uinds, unq_uinv = np.unique(
            self.acoef_final[:, :6], return_index=True,
            return_inverse=True, axis=0,
        )

        # sum the Einstein coefficients for each group of unique upper levels
        # sum_acoefs = np.zeros(unq_uinv.shape[0])
        sum_acoefs = np.zeros(np.unique(unq_uinv).shape)

        for i, ii in enumerate(unq_uinv):
            sum_acoefs[ii] += self.acoef_final[i, -1]

        # the lifetimes are the inverse of the acumulated sums
        lifetimes = 1.0 / sum_acoefs

        # conncatenate the upper levels and the calculated lifetimes
        lifetimes_final = np.column_stack(
            (self.acoef_final[unq_uinds, :6], lifetimes[:unq_uinds.shape[0]])
        )

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

        kt = self.c_boltzmannk * T
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

    def _define_optional_keywords(self):

        return {
            1: 'spectrum_type',
            2: 'apply_rules',
            3: 'compute_A',
            4: 'compute_rel_intensity',
            6: 'temperature',
            7: 'save_hlf',
            8: 'filename_rel_int',
            9: 'filename_Acoefs',
            10: 'filename_hlf',
            11: 'compare_obs',
            12: 'ninter',
            13: 'wavenumber_only',
            14: 'size',
            15: 'filename_wavenumber'
        }

    def compute_linelist(self, Elower=None, Eupper=None,
                         ujrange=None, ljrange=None, lsymm=None, usymm=None,
                         freq_range=None, dmf=None, **kwargs):

        # get optional keywords
        opt_keys = self._define_optional_keywords()

        # spectrum_type = kwargs.get(opt_keys[1]) or 'absorption'
        apply_rules = kwargs.get(opt_keys[2]) or True
        temperature = kwargs.get(opt_keys[6]) or 297
        ninter = kwargs.get(opt_keys[12]) or 1000
        compute_A = kwargs.get(opt_keys[3]) or False
        compute_rel_intensity = kwargs.get(opt_keys[4]) or False
        wavenumber_only = kwargs.get(opt_keys[13]) or False
        save_hlf = kwargs.get(opt_keys[7]) or False
        filename_hlf = kwargs.get(opt_keys[10]) or 'hlf.dat'
        filename_rel_int = kwargs.get(opt_keys[8]) or 'result_hlf.dat'
        filename_Acoefs = kwargs.get(opt_keys[9]) or 'result_A.dat'
        size = kwargs.get(opt_keys[14]) or 10000
        wavenumbers_file = kwargs.get(opt_keys[15]) or 'wavenumbers.dat'

        # TODO: check if dmf is None
        freq_range = freq_range or (0, 1e5)
        usymm = usymm or (0, 1)
        lsymm = lsymm or (0, 1)

        # js = min(np.amin(uterms[:, uj]), np.amin(lterms[:, lj]))
        # je = max(np.amax(uterms[:, uj]), np.amax(lterms[:, lj]))

        evalues = self.mlevels.get_predicted_data()

        # TODO: set default constarints for J, par, El, Eu, freq

        eind, jind, pind = 1, 2, 3

        # constraints by energy
        lower_levels = evalues[
            (evalues[:, eind] >= Elower[0]) & (evalues[:, eind] <= Elower[1])
        ]
        upper_levels = evalues[
            (evalues[:, eind] >= Eupper[0]) & (evalues[:, eind] <= Eupper[1])
        ]

        # constraints by J
        ujs, uje = ujrange[0], ujrange[1]
        # upper_jrots = np.arange(ujs, uje+1)
        upper_levels = upper_levels[
            (upper_levels[:, jind] >= ujs) & (upper_levels[:, jind] <= uje)
        ]

        ljs, lje = ljrange[0], ljrange[1]
        # lower_jrots = np.arange(ljs, lje+1)
        lower_levels = lower_levels[
            (lower_levels[:, jind] >= ljs) & (lower_levels[:, jind] <= lje)
        ]

        # constarints by symmetry
        lower_levels = lower_levels[
            (np.in1d(lower_levels[:, pind], lsymm[0])) |
            (np.in1d(lower_levels[:, pind], lsymm[1]))
        ]
        upper_levels = upper_levels[
            (np.in1d(upper_levels[:, pind], usymm[0])) |
            (np.in1d(upper_levels[:, pind], usymm[1]))
        ]

        freq_calc = self._create_wavenumbers_list(
            upper_levels, lower_levels, jind, pind, freq_range, size
        )

        # extract initial and final levels from the wavenumbers list
        # by filtering the unique inidicies of the levels
        _, ilinds = np.unique(freq_calc[:, 0], return_index=True)
        _, flinds = np.unique(freq_calc[:, 9], return_index=True)
        initial_levels = freq_calc[:, :9]
        final_levels = freq_calc[:, 9:-1]
        print(compute_A, compute_rel_intensity)
        compute_A = compute_A or compute_rel_intensity
        print(compute_A, compute_rel_intensity)
        apply_rules = apply_rules or (compute_A or compute_rel_intensity)

        if wavenumber_only:
            compute_A, compute_rel_intensity, save_hlf = False, False, False
            cols = [6, 2, 1, 3, 4, 5, 15, 11, 10, 12, 13, 14, -1]
            self._save_wavenumbers(freq_calc[:, cols], wavenumbers_file)

        hlf = np.ones(freq_calc.shape[0])
        if apply_rules:
            usymm, lsymm = freq_calc[:, 3], freq_calc[:, 12]
            uj, lj = freq_calc[:, 2], freq_calc[:, 11]
            ulambda, llambda = freq_calc[:, 7], freq_calc[:, -3]
            uomega, lomega = freq_calc[:, 8], freq_calc[:, -2]

            hlf = self._compute_Honl_London_factors(
                usymm, lsymm, uj, lj, uomega, lomega, ulambda, llambda
            )

        if save_hlf:
            calc_hlf = np.c_[freq_calc, hlf]
            # calc_hlf = calc_hlf[np.where(hlf != 0.0)[0], :]
            cols = [6, 2, 1, 3, 4, 5, 15, 11, 10, 12, 13, 14, -2, -1]
            self._save_Honl_London_factor(calc_hlf[:, cols], filename_hlf)

        if compute_A or compute_rel_intensity:

            # calculate Einstein coeffcients
            result = self._compute_electric_dipole_elements(
                initial_levels, final_levels, dmf, ninter=ninter
            )

            line_strength = result[:, -1] * hlf

            jinitial, wavenumber = result[:, 1], result[:, -2]
            gji = (2 * jinitial + 1)  # statistical weight of the initial level
            acoef = self._calculate_Einstein_coeffcients(
                wavenumber, line_strength, gji
            )
            acoef_result = np.c_[result[:, :-1], acoef]
            nonzero_inds = np.where(acoef != 0.0)[0]
            nonzero_result = acoef_result[nonzero_inds, :]

            self._save_Einstein_coeffcients(nonzero_result, filename_Acoefs)

            if compute_rel_intensity:
                kt = self.c_boltzmannk * temperature
                # TODO: check if the exponent should be calculated with
                # the energy value or the term value

                intensity = acoef * gji * np.exp(-freq_calc[:, 10] / kt)
                # (1 - np.exp(-freq_calc[:, -1]/kt))
                # np.square(np.square(freq_calc[:, -1]))
                intensity_result = np.c_[result[:, :-1], intensity]
                nonzero_result = intensity_result[nonzero_inds, :]
                cols = [6, 2, 1, 3, 4, 5, 15, 11, 10, 12, 13, 14, -2, -1]

                self._save_rel_intensity(nonzero_result, filename_rel_int)

    def _create_wavenumbers_list(self, upper_levels, lower_levels, jind, pind,
                                 freq_range, size):

        lshape, ushape = lower_levels.shape, upper_levels.shape
        freq_calc = np.zeros((size, 19))
        countf = 0

        # select which columns from the levels to use
        # n' E' J' p' iso' st' v' l' o' n E J p iso st v l o freq
        select = [0, 1, 2, 3, 4, -4, -3, -2, -1]
        jcol, pcol = jind + len(select), pind + len(select)

        for ulevel in range(ushape[0]):
            level = np.tile(
                upper_levels[ulevel, select], lshape[0]
            ).reshape(lshape[0], len(select))

            levels = np.hstack((level, lower_levels[:, select]))

            # apply selection rules by J and symmetry - the remaining
            # approximate selection rules will be applied by computing HLF
            jp = levels[:, jind] == levels[:, jcol] - 1.0
            jq = levels[:, jind] == levels[:, jcol]
            jr = levels[:, jind] == levels[:, jcol] + 1.0

            cond_pee = jp & (levels[:, pind] == 1) & (levels[:, pcol] == 1)
            cond_pff = jp & (levels[:, pind] == 0) & (levels[:, pcol] == 0)
            cond_qef = jq & (levels[:, pind] == 1) & (levels[:, pcol] == 0)
            cond_qfe = jq & (levels[:, pind] == 0) & (levels[:, pcol] == 1)
            cond_ree = jr & (levels[:, pind] == 1) & (levels[:, pcol] == 1)
            cond_rff = jr & (levels[:, pind] == 0) & (levels[:, pcol] == 0)

            levels_selected = np.vstack((
                levels[cond_pee], levels[cond_pff], levels[cond_qef],
                levels[cond_qfe], levels[cond_ree], levels[cond_rff]
            ))

            freq = levels_selected[:, 1] - levels_selected[:, 10]
            fq_calc = np.column_stack((levels_selected, freq))

            freq_calc[countf:countf+fq_calc.shape[0]] = fq_calc
            countf += fq_calc.shape[0]

        freq_calc = freq_calc[:countf]

        # apply constraints by frequency
        fs, fe = freq_range[0], freq_range[1]
        fcond = (freq_calc[:, -1] >= fs) & (freq_calc[:, -1] <= fe)
        fcond = fcond & (freq_calc[:, -1] >= 0.0)
        freq_calc = freq_calc[fcond]

        if freq_calc.shape[0] == 0:
            raise SystemExit(
                'Error: Allowed transition frequencies were not found.\n'
            )

        return freq_calc

    def _calculate_Einstein_coeffcients(self, wavenumber, line_strength, gji):

        # e0 = constants.value('vacuum electric permittivity')
        # planck = constants.value('Planck constant')
        # consts = (16 * np.pi**3) / 3 * e0 * planck
        # acoef = (consts * wavenumber**3 * line_strength) / gi
        acoef = (3.1361891e-7 * line_strength * wavenumber**3) / gji

        return acoef

    def _compute_electric_dipole_elements(self, init_levels, final_levels,
                                          dmf, ninter=1000):

        ivec_inds = init_levels[:, 0].astype(np.int)
        fvec_inds = final_levels[:, 0].astype(np.int)
        # print('inds', ivec_inds, fvec_inds)

        select = [-3, 2, 1, 3, 4, -4]
        result = np.zeros((ivec_inds.shape[0], 2*len(select)+2))

        rgrid, rstep = self.grid.rgrid * C_bohr, self.grid.rstep * C_bohr
        ngrid = rgrid.shape[0]
        igrid, istep = np.linspace(
            self.grid.rmin, self.grid.rmax, num=ninter,
            endpoint=True, retstep=True
        )
        igrid, istep = igrid * C_bohr, istep * C_bohr

        sinc_matrix = self._calculate_sinc_matrix(rgrid, igrid, ngrid, rstep)
        if dmf is None:
            dipole_matrix = np.ones(
                (self.mlevels.nch, self.mlevels.nch, igrid.shape[0])
            )
        else:
            dipole_matrix = self._interpolate_dipole_moment(rgrid, igrid, dmf)

        # np.savetxt('sinc_matrix.dat', sinc_matrix, fmt='%14.11e')
        # np.savetxt('evecs_all.dat', mlevels.evecs_matrix, fmt='%12.11e')
        # np.savetxt('dipole_matrix.dat', dipole_matrix[0, 0, :], fmt='%14.8e')

        # arrays with the initial and final levels are guaranteed to have the
        # same size since they are extracted from the list of wavenumbers so
        # we can loop over them with only one cycle (no need of 2 nested loops)

        ii = 0
        for (i, j) in zip(ivec_inds, fvec_inds):
            dme = 0.0
            for n in range(self.mlevels.nch):
                for m in range(self.mlevels.nch):
                    ncoefs = self.mlevels.evecs_matrix[n*ngrid:(n+1)*ngrid, i-1]
                    kcoefs = self.mlevels.evecs_matrix[m*ngrid:(m+1)*ngrid, j-1]

                    # for l in range(ngrid):
                    #     for p in range(ngrid):
                    #         sincl = sinc_matrix[l, :]
                    #         sincp = sinc_matrix[p, :]
                    #         sumq = np.sum(sincl*dipole_matrix[n, m, :]*sincp)
                    #         dme += kcoefs[l] * ncoefs[p] * sumq * istep
                    # res = dme / rstep
                    # print(n, m, res**2)
                    # to get the fianl result this should be outside n, m loops
                    # res = dme / rstep

                    sumq = np.dot(
                        sinc_matrix, (dipole_matrix[n, m, :]*sinc_matrix).T
                    )
                    dme += np.dot(kcoefs, ncoefs * sumq) * istep
            res = np.sum(dme) / rstep

            freq = init_levels[ii, 1] - final_levels[ii, 1]
            result[ii, :] = np.hstack((
                init_levels[ii, select], final_levels[ii, select], freq, res**2
            ))
            ii += 1

            # print(i, j, freq, res**2)

        return result

    def _interpolate_dipole_moment(self, rgrid, igrid, dmf):

        # TODO: initilize with ones and remove try/except?
        # TODO: check if getvalue works when dmf is string
        # print(io.StringIO(dmf[(n, k)]))
        # print(io.StringIO(dmf[(n, k)]).getvalue())
        dipole_moment = np.zeros(
            (self.mlevels.nch, self.mlevels.nch, igrid.shape[0])
        )

        for n in range(1, self.mlevels.nch+1):
            for k in range(1, self.mlevels.nch+1):
                try:
                    dip_moment = np.loadtxt(io.StringIO(dmf[(n, k)]).getvalue())
                    dmx, dmy = dip_moment[:, 0], dip_moment[:, 1]
                    # dmx /= C_bohr
                except KeyError:
                    dmx = rgrid
                    dmy = np.ones_like(dmx)

                cs = CubicSpline(dmx, dmy, bc_type='natural')
                dipole_moment[n-1, k-1, :] = cs(igrid)

        return dipole_moment

    def _calculate_sinc_matrix(self, rgrid, igrid, ngrid, rstep):

        sinc_matrix = np.zeros((ngrid, igrid.shape[0]))

        for i in range(ngrid):
            argi = (igrid - rgrid[i]) / rstep
            sinc = np.sinc(argi)
            sinc_matrix[i, :] = sinc

        return sinc_matrix

    def _compute_Honl_London_factors(self, usymm, lsymm, uj, lj, uomega,
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

    def _save_wavenumbers(self, freq_calc, filename):

        labels = (
            "v'", "J'", "E'", "symm'", "iso'", "state'", 'v',
            'J', 'E', 'symm', 'iso', 'state', 'freq'
        )
        header = (
            f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<6}'
            f'{labels[4]:<6}{labels[5]:<9}{labels[6]:<6}{labels[7]:<13}'
            f'{labels[8]:<9}{labels[9]:<6}{labels[10]:<6}{labels[11]:<10}'
            f'{labels[12]}'
        )
        fmt = (
            '%6.1d', '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%7.1d',
            '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%15.5f'
        )
        np.savetxt(filename, freq_calc, header=header, fmt=fmt)

    def _save_Honl_London_factor(self, hlf, filename):

        labels = (
            "v'", "J'", "E'", "symm'", "iso'", "state'", 'v',
            'J', 'E', 'symm', 'iso', 'state', 'freq', 'hlf'
        )
        header = (
            f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<6}'
            f'{labels[4]:<6}{labels[5]:<9}{labels[6]:<6}{labels[7]:<13}'
            f'{labels[8]:<9}{labels[9]:<6}{labels[10]:<6}{labels[11]:<10}'
            f'{labels[12]:<14}{labels[13]}'
        )
        fmt = (
            '%6.1d', '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%7.1d',
            '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%15.5f', '%12.5f'
        )
        np.savetxt(filename, hlf, header=header, fmt=fmt)

    def _save_rel_intensity(self, result, out_file):

        labels = (
            "v'", "J'", "E'", "symm'", "iso'", "state'", 'v',
            'J', 'E', 'symm', 'iso', 'state', 'freq', 'intensity'
        )
        header = (
            f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<6}'
            f'{labels[4]:<6}{labels[5]:<9}{labels[6]:<6}{labels[7]:<13}'
            f'{labels[8]:<9}{labels[9]:<6}{labels[10]:<6}{labels[11]:<10}'
            f'{labels[12]:<15}{labels[13]}'
        )
        fmt = (
            '%5.1d', '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%7.1d',
            '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%15.5f', '%16.5e',
        )
        np.savetxt(out_file, result, comments='#', header=header, fmt=fmt)

    def _save_Einstein_coeffcients(self, acoefs, filename):

        labels = (
            "v'", "J'", "E'", "symm'", "iso'", "state'", 'v',
            'J', 'E', 'symm', 'iso', 'state', 'freq', 'A'
        )
        header = (
            f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<6}'
            f'{labels[4]:<6}{labels[5]:<9}{labels[6]:<6}{labels[7]:<13}'
            f'{labels[8]:<9}{labels[9]:<6}{labels[10]:<6}{labels[11]:<10}'
            f'{labels[12]:<17}{labels[13]}'
        )
        fmt = (
            '%5.1d', '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%7.1d',
            '%7.1f', '%14.5f', '%5.1d', '%5.1d', '%5.1d', '%15.5f', '%15.5e'
        )
        np.savetxt(filename, acoefs, comments='#', header=header, fmt=fmt)

    def _save_lifetimes(self, lifetimes, filename):

        labels = (
            "v'", "J'", "E'", "symm'", "iso'", "state'", 'life_time'
        )
        header = (
            f'{labels[0]:^9}{labels[1]:<13}{labels[2]:<9}{labels[3]:<8}'
            f'{labels[4]:<7}{labels[5]:<10}{labels[6]}'
        )
        fmt = (
            '%5.1d', '%7.1f', '%14.5f', '%6.1d', '%6.1d', '%6.1d', '%17.5e'
        )
        np.savetxt(filename, lifetimes, comments='#', header=header, fmt=fmt)

    def _get_default_jrange(self, uterms, uj, lterms, lj):

        js = min(np.amin(uterms[:, uj]), np.amin(lterms[:, lj]))
        je = max(np.amax(uterms[:, uj]), np.amax(lterms[:, lj]))
        return js, je

    def _save_compared_frequencies(self, fname, final_freqs, fmode):

        labels = (
            "v'", "v", "J'", "J", "st'", "st", "p'", "p", "FJ'",
            "FJ", "calc_freq", "obs_freq", "calc-obs", "unc",
        )

        header = (
            f'{labels[0]:^10}{labels[1]:<7}{labels[2]:<8}{labels[3]:<6}'
            f'{labels[4]:<6}{labels[5]:<7}{labels[6]:<6}{labels[7]:<10}'
            f'{labels[8]:<12}{labels[9]:<10}{labels[10]:<12}{labels[11]:<13}'
            f'{labels[12]:<11}{labels[13]}'
        )

        if fmode == 'ab':
            header = ''

        fmt = [
            '%6d', '%5d', '%7.1f', '%7.1f', '%5d', '%5d', '%5d', '%5d',
            '%13.4f', '%12.4f', '%12.4f', '%12.4f', '%9.4f', '%8.3f'
        ]

        with open(fname, fmode) as f:
            np.savetxt(f, final_freqs, header=header, comments='#', fmt=fmt)

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

    def _compare_frequencies(self, freq_exp_file, fcomp_file, freq_calc):

        assigned_freq = np.loadtxt(freq_exp_file)
        fmode = 'w'

        for row in range(0, freq_calc.shape[0]):
            ln = assigned_freq.shape[0]
            term = np.tile(
                freq_calc[row, :], ln
            ).reshape(ln, freq_calc.shape[1])
            merged_freq = np.hstack((term, assigned_freq))

            # if freq_calc.shape[0] > assigned_freq.shape[0]:
            #     merged_freq = np.hstack(
            #         (freq_calc[0:assigned_freq.shape[0], :], assigned_freq)
            #     )
            # else:
            #     merged_freq = np.hstack(
            #         (freq_calc, assigned_freq[0:freq_calc.shape[0], :])
            #     )

            vcond = (merged_freq[:, 0] == merged_freq[:, 16]) & \
                    (merged_freq[:, 1] == merged_freq[:, 17])
            merged_freq = merged_freq[vcond]

            scond = (merged_freq[:, 8] == merged_freq[:, 24]) & \
                    (merged_freq[:, 9] == merged_freq[:, 25])
            merged_freq = merged_freq[scond]

            rcond = (merged_freq[:, 2] == merged_freq[:, 18]) & \
                    (merged_freq[:, 3] == merged_freq[:, 19])
            merged_freq = merged_freq[rcond]

            pcond = (merged_freq[:, 10] == merged_freq[:, 26]) & \
                    (merged_freq[:, 11] == merged_freq[:, 27])
            merged_freq = merged_freq[pcond]

            # lcond = (merged_freq[:, 4] == merged_freq[:, 20]) & \
            #         (merged_freq[:, 5] == merged_freq[:, 21])
            # merged_freq = merged_freq[lcond]

            # ocond = (merged_freq[:, 6] == merged_freq[:, 22]) & \
            #         (merged_freq[:, 7] == merged_freq[:, 23])
            # merged_freq = merged_freq[ocond]

            # will work only for the if case above
            merged_freq_final = np.column_stack((
                merged_freq[:, 0:4], merged_freq[:, 8:12],
                merged_freq[:, 12:15], merged_freq[:, 28],
                merged_freq[:, 14] - merged_freq[:, 28],
                merged_freq[:, 29]
            ))

            if merged_freq_final.shape[0] != 0:
                #     branches = \
                #         np.vectorize(self._map_number_to_branch().get)(freqs[:, -1])
                #     branch_freqs = np.column_stack((freqs[:, :-1], branches))
                # TODO: add branch inside this function
                # self._save_compared_frequencies(fcomp_file, branch_freqs)
                self._save_compared_frequencies(fcomp_file, merged_freq_final, fmode)
                fmode = 'ab'
            # else:
            #     print(
            #         f'No calculated frequencies corresponding to the '
            #         f'observed frequencies in file - {freq_exp_file} found.'
            #     )
