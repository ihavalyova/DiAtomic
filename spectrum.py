import os
import math
import py3nj
import numpy as np
# from scipy.interpolate import CubicSpline


class Spectrum:
    def __init__(self, dmf, spec_type='absorption', wavefunc1=None):

        self.spec_type = spec_type

        self.calc_freq_file = 'calculated_frequencies.dat'
        self.comp_freq_file = 'compared_frequencies.dat'
        self.evals_file = 'eigenvalues_all.dat'

        if self.spec_type not in ('absorption', 'emission', 'a', 'e'):
            raise SystemExit(
                f'"{spec_type}" is not a valid spectrum type'
            )

    def calculate_frequencies(self, **kwargs):

        """Will calculate the possible transition frequencies
        for one of the two options:
        (i) between two states
        (ii) between arbitrary number of prevously computed term values

        Raises:
            SystemExit: if the channel does not match the required type
            SystemExit: if the files with term values are not found
            SystemExit: if the channels or terms are not provided in pairs

        Remarks:
            0. by channels and by terms are two separate options and cannot
               be used simultaneously
            1. uch/lch should be integer numbers
            2. uterms/lterms are names of an existing files in specific format
            3. if channels and terms files are specified simultaneously,
               the default channel option will be used
            4. if uch is specified lch should also be specified;
               the same applies for lterms and uterms
        """

        uch = kwargs.get('uch')
        lch = kwargs.get('lch')
        uterms = kwargs.get('uterms')
        lterms = kwargs.get('lterms')

        if uch is not None and lch is not None:
            try:
                self.check_channels_type(uch)
                self.check_channels_type(lch)
            except TypeError as te:
                raise SystemExit(te)

            self.calculate_frequencies_by_channels(uch, lch, kwargs)
        elif uterms is not None and lterms is not None:
            try:
                self.check_terms_file_type(uterms)
                self.check_terms_file_type(lterms)
            except FileNotFoundError as fe:
                raise SystemExit(fe)

            self.calculate_frequencies_by_term_values(uterms, lterms, kwargs)
        else:
            raise SystemExit(
                f'Invalid input in Spectrum object: '
                f'lch={lch}, uch={uch} and uterms={uterms}, lterms={lterms} '
                f'should be specified in pairs but not simultaneously.'
            )

    def check_channels_type(self, ch):

        if not isinstance(ch, (int)):
            raise TypeError(
                f'Invalid type in spectrum: {ch} should be an integer number'
            )

    def check_terms_file_type(self, terms):

        if not os.path.isfile(terms):
            raise FileNotFoundError(
                f'Terms={terms} should be a name of an existing file.'
            )

    def define_keywords(self):

        return {
            0: 'evalues_file',
            1: 'output_freq_file',
            2: 'obs_freq_file',
            3: 'compared_freq_file',
            4: 'Jrange',
            5: 'vrange',
            6: 'uparity',
            7: 'lparity',
            8: 'upper_state',
            9: 'lower_state',
            10: 'freq_range',
            11: 'HLF',
            12: 'extract_nonzero',
            13: 'grid',
            14: 'isotopes',
            15: 'FCF',
            16: 'bands',
        }

    def calculate_frequencies_by_channels(self, uch, lch, kwargs):

        """ will compute all allowed transitions between two channels
        for which the eigenvalues and eigenvectors have been calculated

        Args:
            uch (int): the number of the upper channel/state
            lch (int): the number of the lower channel/state
            kwargs (dict): the provided options as keywords

        Remarks:
            1. The keywords 'upper_state' and 'lower_state' are
            not relevant in this case. They are used in the computation
            by term values when more then two states can be invloved
            2. The keyword 'obs_freq_file' should be a name of an existing file
            3. When FCF should be computed grid option is required
            4. When keyword 'HLF'/'FCF' is not set or set to False,
            the Honl-London/Frank-Condon factors will not be computed
            and written in the output files
            5. The type of bands is tuple of tuples or list of lists
        """

        keys = self.define_keywords()

        evalues_file = kwargs.get(keys[0]) or self.evals_file
        evalues = np.loadtxt(evalues_file)

        uterms = evalues[evalues[:, -3] == uch]
        uterms = uterms[:, [0, 1, 2, 3, 4, -3, -2, -1]]
        lterms = evalues[evalues[:, -3] == lch]
        lterms = lterms[:, [0, 1, 2, 3, 4, -3, -2, -1]]

        predicted_freq_file = kwargs.get(keys[1]) or self.calc_freq_file
        obs_freq_file = kwargs.get(keys[2])
        compared_freq_file = kwargs.get(keys[3]) or self.comp_freq_file

        jrange = kwargs.get(keys[4])
        vrange = kwargs.get(keys[5])
        upar = kwargs.get(keys[6]) or (0, 1)
        lpar = kwargs.get(keys[7]) or (0, 1)
        frange = kwargs.get(keys[10]) or (-np.inf, np.inf)
        compute_hlf = kwargs.get(keys[11]) or False
        extract = kwargs.get(keys[12]) or False

        # TODO: check id grid is provided when computing FCF
        grid = kwargs.get(keys[13])

        isotopes = kwargs.get(keys[14]) or [1]
        compute_fcf = kwargs.get(keys[15]) or False
        bands = kwargs.get(keys[16]) or ((0, 0),)

        js, je = self.get_jrange_values(jrange, uterms, lterms)
        vs, ve = self.get_vrange_values(vrange, uterms, lterms)
        fs, fe = frange[0], frange[1]

        freq_calc = np.array([])

        for row in range(0, uterms.shape[0]):

            uterm = np.tile(uterms[row, :], lterms.shape[0]). \
                reshape(lterms.shape[0], uterms.shape[1])

            d = np.hstack((uterm, lterms))

            vcond = \
                (d[:, 1] >= vs) & (d[:, 1] <= ve) & \
                (d[:, 9] >= vs) & (d[:, 9] <= ve)
            d = d[vcond]

            jcond = \
                (d[:, 3] >= js) & (d[:, 3] <= je) & \
                (d[:, 11] >= js) & (d[:, 11] <= je)
            d = d[jcond]

            pcond = (np.in1d(d[:, 4], upar)) & (np.in1d(d[:, 12], lpar))
            d = d[pcond]

            enr_cond = (d[:, 2] >= d[:, 10])
            jp = d[:, 3] == d[:, 11] - 1.0
            jq = d[:, 3] == d[:, 11]
            jr = d[:, 3] == d[:, 11] + 1.0

            cond_pee = jp & (d[:, 4] == 1) & (d[:, 12] == 1) & enr_cond
            cond_pff = jp & (d[:, 4] == 0) & (d[:, 12] == 0) & enr_cond
            cond_qef = jq & (d[:, 4] == 1) & (d[:, 12] == 0) & enr_cond
            cond_qfe = jq & (d[:, 4] == 0) & (d[:, 12] == 1) & enr_cond
            cond_ree = jr & (d[:, 4] == 1) & (d[:, 12] == 1) & enr_cond
            cond_rff = jr & (d[:, 4] == 0) & (d[:, 12] == 0) & enr_cond

            data_pee = d[cond_pee]
            data_pff = d[cond_pff]
            data_qef = d[cond_qef]
            data_qfe = d[cond_qfe]
            data_ree = d[cond_ree]
            data_rff = d[cond_rff]

            shapes = np.array([
                data_pee.shape[0], data_pff.shape[0],
                data_qef.shape[0], data_qfe.shape[0],
                data_ree.shape[0], data_rff.shape[0]
            ])

            # allowed frequencies

            fq_pee = (data_pee[:, 2] - data_pee[:, 10])[:, np.newaxis]
            fq_pff = (data_pff[:, 2] - data_pff[:, 10])[:, np.newaxis]
            fq_qef = (data_qef[:, 2] - data_qef[:, 10])[:, np.newaxis]
            fq_qfe = (data_qfe[:, 2] - data_qfe[:, 10])[:, np.newaxis]
            fq_ree = (data_ree[:, 2] - data_ree[:, 10])[:, np.newaxis]
            fq_rff = (data_rff[:, 2] - data_rff[:, 10])[:, np.newaxis]

            # branches

            brnach_pee = np.repeat(np.array([1]), shapes[0])[:, np.newaxis]
            brnach_pff = np.repeat(np.array([2]), shapes[1])[:, np.newaxis]
            brnach_qef = np.repeat(np.array([3]), shapes[2])[:, np.newaxis]
            brnach_qfe = np.repeat(np.array([4]), shapes[3])[:, np.newaxis]
            brnach_ree = np.repeat(np.array([5]), shapes[4])[:, np.newaxis]
            brnach_rff = np.repeat(np.array([6]), shapes[5])[:, np.newaxis]

            data_pee = np.hstack((data_pee, np.hstack((fq_pee, brnach_pee))))
            data_pff = np.hstack((data_pff, np.hstack((fq_pff, brnach_pff))))
            data_qef = np.hstack((data_qef, np.hstack((fq_qef, brnach_qef))))
            data_qfe = np.hstack((data_qfe, np.hstack((fq_qfe, brnach_qfe))))
            data_ree = np.hstack((data_ree, np.hstack((fq_ree, brnach_ree))))
            data_rff = np.hstack((data_rff, np.hstack((fq_rff, brnach_rff))))

            freq_calc = np.append(freq_calc, np.vstack((
                data_pee, data_pff, data_qef,
                data_qfe, data_ree, data_rff
            )))

        freq_calc = freq_calc.reshape(
            -1, uterms.shape[1] + lterms.shape[1] + 2
        )

        s = freq_calc.shape[1]
        freq_calc = freq_calc[
            :, [1, 9, 3, 11, 6, 14, 7, 15, 5, 13, 4, 12, 2, 10, s-2, s-1]
        ]

        fcond = (freq_calc[:, 14] >= fs) & (freq_calc[:, 14] <= fe)
        freq_calc = freq_calc[fcond]

        branch_freq = self.find_branch_frequencies(freq_calc, compute_hlf)

        self.save_computed_frequencies(
            predicted_freq_file, branch_freq, compute_hlf
        )

        if extract and compute_hlf:
            hlf = branch_freq[:, -1].astype(float)

            # find the indices of the elements with nonzero Honl-London factors
            inonzero = np.where(hlf != 0.0)[0]

            extract_file = 'extract_' + predicted_freq_file

            self.save_computed_frequencies(
                extract_file, branch_freq[inonzero], compute_hlf
            )

        if obs_freq_file is not None:
            self.compare_frequency_lists(
                obs_freq_file, compared_freq_file, freq_calc
            )

        if compute_fcf:
            fcfs = self.compute_Frank_Condon_factors(
                uch, lch, grid, (js, je), isotopes, bands
            )

            fcf_file = 'FCF_' + predicted_freq_file

            self.save_computed_FCFs(fcf_file, fcfs)

    def find_branch_frequencies(self, freq_calc, compute_hlf):

        branches = \
            np.vectorize(self.map_number_to_branch().get)(freq_calc[:, -1])
        branch_freq = np.column_stack((freq_calc[:, :-1], branches))

        hlf = None
        if compute_hlf:
            hlf = self.compute_Honl_London_factors(freq_calc)
            branch_freq = np.column_stack((branch_freq, hlf))

        return branch_freq

    def compare_frequency_lists(self, obs_freq_file, compared_freq_file,
                                freq_calc):

        freqs = self.compare_frequencies(obs_freq_file, freq_calc)

        if freqs.shape[0] != 0:
            branches = \
                np.vectorize(self.map_number_to_branch().get)(freqs[:, -1])

            branch_freqs = np.column_stack((freqs[:, :-1], branches))

            self.save_compared_frequencies(compared_freq_file, branch_freqs)
        else:
            print(
                f'No calculated frequencies corresponding to the '
                f'observed frequencies in file - {obs_freq_file} found.'
            )

    def get_jrange_values(self, jrange, uterms, lterms):

        if jrange is None:
            js = min(np.amin(uterms[:, 2]), np.amin(lterms[:, 2]))
            je = max(np.amax(uterms[:, 2]), np.amax(lterms[:, 2]))
            return js, je

        return jrange[0], jrange[1]

    def get_vrange_values(self, vrange, uterms, lterms):

        if vrange is None:
            vs = min(np.amin(uterms[:, 1]), np.amin(lterms[:, 1]))
            ve = max(np.amax(uterms[:, 1]), np.amax(lterms[:, 1]))
            return vs, ve

        return vrange[0], vrange[1]

    # TODO: change the indices in this function since the
    # format of eigenvalues_all.dat file was changed
    def calculate_frequencies_by_term_values(self, uterms, lterms, kwargs):

        keys = self.define_keywords()

        terms = np.loadtxt(uterms)  # upper terms
        evalues = np.loadtxt(lterms)  # lower terms

        predicted_freq_file = kwargs.get(keys[1]) or self.calc_freq_file

        # 'obs_freq_file' should be a name of an existing file
        obs_freq_file = kwargs.get(keys[2])
        compared_freq_file = kwargs.get(keys[3]) or self.comp_freq_file

        jrange = kwargs.get(keys[4])
        vrange = kwargs.get(keys[5])
        upar = kwargs.get(keys[6])
        lpar = kwargs.get(keys[7])
        su = np.array(kwargs.get(keys[8]))
        sl = np.array(kwargs.get(keys[9]))
        frange = kwargs.get(keys[10])
        compute_hlf = kwargs.get(keys[11]) or False
        extract = kwargs.get(keys[12]) or False

        js, je = self.get_jrange_values(jrange, terms, evalues)
        vs, ve = self.get_vrange_values(vrange, terms, evalues)
        fs, fe = self.get_frange_values(frange)

        freq_calc = np.array([])

        for row in range(0, terms.shape[0]):
            term = np.tile(terms[row, :], evalues.shape[0]). \
                reshape(evalues.shape[0], terms.shape[1])

            d = np.hstack((term, evalues))

            vcond = \
                (d[:, 1] >= vs) & (d[:, 1] <= ve) & \
                (d[:, 9] >= vs) & (d[:, 9] <= ve)
            jcond = \
                (d[:, 2] >= js) & (d[:, 2] <= je) & \
                (d[:, 10] >= js) & (d[:, 10] <= je)
            scond = (np.in1d(d[:, 7], su)) & (np.in1d(d[:, 15], sl))
            pcond = (np.in1d(d[:, 6], upar)) & (np.in1d(d[:, 14], lpar))

            d = d[vcond]
            d = d[jcond]
            d = d[scond]
            d = d[pcond]

            enr_cond = (d[:, 5] >= d[:, 13])
            jp = d[:, 2] == d[:, 10] - 1.0
            jq = d[:, 2] == d[:, 10]
            jr = d[:, 2] == d[:, 10] + 1.0

            cond_pee = jp & (d[:, 6] == 1) & (d[:, 14] == 1) & enr_cond
            cond_pff = jp & (d[:, 6] == 0) & (d[:, 14] == 0) & enr_cond
            cond_qef = jq & (d[:, 6] == 1) & (d[:, 14] == 0) & enr_cond
            cond_qfe = jq & (d[:, 6] == 0) & (d[:, 14] == 1) & enr_cond
            cond_ree = jr & (d[:, 6] == 1) & (d[:, 14] == 1) & enr_cond
            cond_rff = jr & (d[:, 6] == 0) & (d[:, 14] == 0) & enr_cond

            data_pee = d[cond_pee]
            data_pff = d[cond_pff]
            data_qef = d[cond_qef]
            data_qfe = d[cond_qfe]
            data_ree = d[cond_ree]
            data_rff = d[cond_rff]

            shapes = np.array([
                data_pee.shape[0], data_pff.shape[0],
                data_qef.shape[0], data_qfe.shape[0],
                data_ree.shape[0], data_rff.shape[0]
            ])

            # allowed frequencies

            fq_pee = (data_pee[:, 5] - data_pee[:, 13])[:, np.newaxis]
            fq_pff = (data_pff[:, 5] - data_pff[:, 13])[:, np.newaxis]
            fq_qef = (data_qef[:, 5] - data_qef[:, 13])[:, np.newaxis]
            fq_qfe = (data_qfe[:, 5] - data_qfe[:, 13])[:, np.newaxis]
            fq_ree = (data_ree[:, 5] - data_ree[:, 13])[:, np.newaxis]
            fq_rff = (data_rff[:, 5] - data_rff[:, 13])[:, np.newaxis]

            # branches

            brnach_pee = np.repeat(np.array([1]), shapes[0])[:, np.newaxis]
            brnach_pff = np.repeat(np.array([2]), shapes[1])[:, np.newaxis]
            brnach_qef = np.repeat(np.array([3]), shapes[2])[:, np.newaxis]
            brnach_qfe = np.repeat(np.array([4]), shapes[3])[:, np.newaxis]
            brnach_ree = np.repeat(np.array([5]), shapes[4])[:, np.newaxis]
            brnach_rff = np.repeat(np.array([6]), shapes[5])[:, np.newaxis]

            data_pee = np.hstack((data_pee, np.hstack((fq_pee, brnach_pee))))
            data_pff = np.hstack((data_pff, np.hstack((fq_pff, brnach_pff))))
            data_qef = np.hstack((data_qef, np.hstack((fq_qef, brnach_qef))))
            data_qfe = np.hstack((data_qfe, np.hstack((fq_qfe, brnach_qfe))))
            data_ree = np.hstack((data_ree, np.hstack((fq_ree, brnach_ree))))
            data_rff = np.hstack((data_rff, np.hstack((fq_rff, brnach_rff))))

            freq_calc = np.append(freq_calc, np.vstack((
                data_pee, data_pff, data_qef,
                data_qfe, data_ree, data_rff
            )))

        freq_calc = freq_calc.reshape(
            -1, terms.shape[1] + evalues.shape[1] + 2
        )
        s = freq_calc.shape[1]
        freq_calc = freq_calc[
            :, [1, 9, 2, 10, 3, 11, 4, 12, 7, 15, 6, 14, 5, 13, s-2, s-1]
        ]

        fcond = (freq_calc[:, 14] >= fs) & (freq_calc[:, 14] <= fe)
        freq_calc = freq_calc[fcond]

        branch_freq = self.find_branch_frequencies(freq_calc, compute_hlf)

        self.save_computed_frequencies(
            predicted_freq_file, branch_freq, compute_hlf
        )

        if extract and compute_hlf:
            hlf = branch_freq[:, -1].astype(float)

            # find the indices of the elements with nonzero Honl-London factors
            inonzero = np.where(hlf != 0.0)[0]

            extract_file = 'extract_' + predicted_freq_file

            self.save_computed_frequencies(
                extract_file, branch_freq[inonzero], compute_hlf
            )

        if obs_freq_file is not None:
            self.compare_frequency_lists(
                obs_freq_file, compared_freq_file, freq_calc
            )

    def save_computed_frequencies(self, fname, branch_freq, compute_hlf):

        names = [
            "v'", "v", "J'", "J", "Lambda'", "Lambda",
            "Omega'", "Omega", "state'", "state", "p'",
            "p", "FJ'", "FJ", "freq", "branch"
        ]

        fmt = [
            '%3.1s', '%6.1s', '%7.5s', '%7.5s', '%7.1s', '%7.1s',
            '%7.3s', '%7.3s', '%7.1s', '%7.1s', '%7.1s', '%7.1s',
            '%13.10s', '%13.10s', '%13.10s', '%6s'
        ]

        if compute_hlf:
            names.append('HLF')
            fmt.append('%11.8s')

        header = \
            ''.join(['{:>3}'.format(k) for k in names[:1]]) + \
            ''.join(['{:>6}'.format(k) for k in names[1:2]]) + \
            ''.join(['{:>8}'.format(k) for k in names[2:4]]) + \
            ''.join(['{:>9}'.format(k) for k in names[4:8]]) + \
            ''.join(['{:>7}'.format(k) for k in names[8:12]]) + \
            ''.join(['{:>12}'.format(k) for k in names[12:15]]) + \
            ''.join(['{:>14}'.format(k) for k in names[15:16]]) + \
            ''.join(['{:>8}'.format(k) for k in names[16:]])

        np.savetxt(fname, branch_freq, header=header, comments='#', fmt=fmt)

    def save_compared_frequencies(self, fname, branch_freqs):

        names = \
            "v'", "v", "J'", "J", "Lam'", "Lam", \
            "Omega'", "Omega", "state'", "state", \
            "p'", "p", "FJ'", "FJ", " freq_calc", \
            "freq_obs", "freq_delta", "unc", "branch"

        header = \
            ''.join(['{:>6}'.format(k) for k in names[:2]]) + \
            ''.join(['{:>8}'.format(k) for k in names[2:7]]) + \
            ''.join(['{:>8}'.format(k) for k in names[7:12]]) + \
            ''.join(['{:>12}'.format(k) for k in names[12:14]]) + \
            ''.join(['{:>15}'.format(k) for k in names[14:17]]) + \
            ''.join(['{:>9}'.format(k) for k in names[17:]])

        fmt = [
            '%7.1s', '%6.1s', '%7.3s', '%7.3s', '%7.1s', '%7.1s',
            '%7.3s', '%7.3s', '%7.1s', '%7.1s', '%6.1s',
            '%7.1s', '%13.10s', '%13.10s', '%13.10s',
            '%13.10s', '%13.10s', '%10.7s', '%6s'
        ]

        np.savetxt(
            fname,
            branch_freqs,
            header=header,
            comments='#',
            fmt=fmt
        )

    def map_number_to_branch(self):

        return {
            1: 'Pee',
            2: 'Pff',
            3: 'Qef',
            4: 'Qfe',
            5: 'Ree',
            6: 'Rff'
        }

    def map_branch_to_number(self):

        return {
            'Pee': 1,
            'Pff': 2,
            'Qef': 3,
            'Qfe': 4,
            'Ree': 5,
            'Rff': 6
        }

    def compare_frequencies(self, freq_exp_file, freq_calc):

        assigned_freq = np.loadtxt(freq_exp_file)

        assigned_freq = assigned_freq[np.lexsort((
            assigned_freq[:, 12], assigned_freq[:, 11],
            assigned_freq[:, 10], assigned_freq[:, 9],
            assigned_freq[:, 8], assigned_freq[:, 7],
            assigned_freq[:, 6], assigned_freq[:, 5],
            assigned_freq[:, 4], assigned_freq[:, 3],
            assigned_freq[:, 2], assigned_freq[:, 1],
            assigned_freq[:, 0]
        ))]

        freq_calc = freq_calc[np.lexsort((
            freq_calc[:, 14], freq_calc[:, 11], freq_calc[:, 10],
            freq_calc[:, 9], freq_calc[:, 8], freq_calc[:, 7],
            freq_calc[:, 6], freq_calc[:, 5], freq_calc[:, 4],
            freq_calc[:, 3], freq_calc[:, 2], freq_calc[:, 1],
            freq_calc[:, 0]
        ))]

        if freq_calc.shape[0] > assigned_freq.shape[0]:
            merged_freq = np.hstack(
                (freq_calc[0:assigned_freq.shape[0], :], assigned_freq)
            )
        else:
            merged_freq = np.hstack(
                (freq_calc, assigned_freq[0:freq_calc.shape[0], :])
            )

        scond = (merged_freq[:, 8] == merged_freq[:, 24]) & \
                (merged_freq[:, 9] == merged_freq[:, 25])
        rcond = (merged_freq[:, 2] == merged_freq[:, 18]) & \
                (merged_freq[:, 3] == merged_freq[:, 19])
        pcond = (merged_freq[:, 10] == merged_freq[:, 26]) & \
                (merged_freq[:, 11] == merged_freq[:, 27])
        vcond = (merged_freq[:, 0] == merged_freq[:, 16]) & \
                (merged_freq[:, 1] == merged_freq[:, 17])
        lcond = (merged_freq[:, 4] == merged_freq[:, 20]) & \
                (merged_freq[:, 5] == merged_freq[:, 21])
        ocond = (merged_freq[:, 6] == merged_freq[:, 22]) & \
                (merged_freq[:, 7] == merged_freq[:, 23])

        conditions = scond & rcond & pcond & vcond & lcond & ocond

        merged_freq = merged_freq[conditions]

        merged_freq = np.column_stack((
            merged_freq[:, 0:15],
            merged_freq[:, 28],
            merged_freq[:, 14] - merged_freq[:, 28],
            merged_freq[:, 29],
            merged_freq[:, 15]
        ))

        return merged_freq

    def compute_Honl_London_factors(self, freqc):

        ueps = self.map_parity_to_epsilon(freqc[:, 10])
        leps = self.map_parity_to_epsilon(freqc[:, 11])

        uj, lj = freqc[:, 2], freqc[:, 3]
        uomega, lomega = freqc[:, 6], freqc[:, 7]
        ulambda, llambda = freqc[:, 4], freqc[:, 5]

        eps_expr = 1.0 + (ueps * leps * np.power(-1, 1.0 + uj - lj))

        udelta = self.map_quantum_number_to_delta(ulambda)
        udelta *= self.map_quantum_number_to_delta(uomega)

        ldelta = self.map_quantum_number_to_delta(llambda)
        ldelta *= self.map_quantum_number_to_delta(lomega)

        delta_expr = 1.0 + udelta + ldelta - 2.0 * udelta * ldelta
        j_expr = ((2.0 * uj) + 1.0) * ((2.0 * lj + 1.0))

        two_l1 = np.int64(2*uj)
        two_l2 = 2*np.ones(freqc.shape[0], dtype=np.int64)
        two_l3 = np.int64(2*lj)
        two_m1 = np.int64(-2*uomega)
        two_m2 = np.int64(2*(ulambda-llambda))
        two_m3 = np.int64(2*lomega)

        # allow for ambiguous sign in the 3j symbol when one of the Lambda
        # quantum numbers is zero (Ref: J.K.G. Watson, JMS 252 (2008))

        if (ulambda.any() == 0 and llambda.any() != 0) or \
           (ulambda.any() != 0 and llambda.any() == 0):
            two_m2 = np.int64(2*(ulambda+llambda))
            two_m3 = np.int64(-2*lomega)

        # set to zero the qunatum numbers that
        # do not satisfy the following conditions

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

    def map_parity_to_epsilon(self, pars):

        return np.array(list(map(lambda z: -1 if z == 0 else z, pars)))

    def map_quantum_number_to_delta(self, qlist):

        return np.array(list(map(lambda z: 1 if z == 0 else 0, qlist)))

    def compute_Frank_Condon_factors(self, uch, lch, grid,
                                     jrange, isotopes, bands):

        """Will calculate the Frank-Condon factors for P, Q and
        R branches for one or several bands if the eigenvectors
        are previously computed and stored in files. Inside these
        files the eigenvectors for all calculated channels are written.
        They are arranged in columns and each column is associated
        and labeled with given vibrational quantum number and state
        determined by the mixing coefficints

        Args:
            uch (int): upper channel's number
            lch (int): lower channel's number
            grid (obj): grid object
            jrange (iterable): the range of J values to be computed
            isotopes (int): the isotope's number
            bands (iterable): the list with bands to be computed

        Raises:
            SystemExit: when no eigenvector files are found

        Returns:
            numpy ndarray: the computed Frank-Condon factors for
            each branch and band

        Remarks:
            1. always will calculate for both parities of upper and lower state
            2. default band is (0, 0)
            3. default number of interp. points for the wavefunctions is 200
        """

        pars = (0, 1)
        js = math.floor(jrange[0])
        je = math.floor(jrange[1])

        rgrid = grid.rgrid
        shape = rgrid.shape[0]

        ninter = 200
        igrid, istep = np.linspace(
            rgrid[0],
            rgrid[shape-1],
            num=ninter,
            endpoint=True,
            retstep=True
        )

        evec_dir = 'eigenvectors'

        # the names of all eigenvector files
        evec_files_all = os.listdir(evec_dir)

        for iso in isotopes:
            # the names of eigenvector files to be used in the computation
            evec_files = [
                f'evector_J{str(j)}_p{str(par)}_i{str(iso)}.dat'
                for par in pars for j in range(js, je)
            ]

            if len(list(set(evec_files).intersection(evec_files_all))) == 0:
                raise SystemExit(
                    f'Cannot find the required eigenvector files in {evec_dir}'
                )

            fcfs = np.array([])
            counter = 0

            for band in bands:
                uv, lv = band

                for evec_file in evec_files:

                    path1 = os.path.join(evec_dir, evec_file)

                    ''' Qef/Qfe branch '''

                    j1 = int(evec_file.rsplit('_')[1][1:])
                    j2 = j1

                    p1 = int(evec_file.rsplit('_')[2][1:])
                    p2 = 1 - p1

                    # lower state
                    evec_file2 = evec_file.replace('p'+str(p1), 'p'+str(p2))
                    path2 = os.path.join(evec_dir, evec_file2)
                    fcfq = -1

                    if os.path.isfile(path1) and os.path.isfile(path2):
                        fcfq = self.compute_single_FCF(
                            path1, path2, uch, lch, uv, lv,
                            shape, rgrid, igrid, istep
                        )

                        branch = self.map_parity_to_branch('Q', p1, p2)
                        qbranch12 = np.array([
                            uv, lv, j1, j2, uch, lch, p1, p2, iso, fcfq, branch
                        ])
                        fcfs = np.hstack((fcfs, qbranch12))
                        counter += 1

                    ''' Rff/Ree branch '''

                    j1 = int(evec_file.rsplit('_')[1][1:])
                    j2 = j1 - 1

                    p1 = int(evec_file.rsplit('_')[2][1:])
                    p2 = p1

                    # lower wavefunction
                    evec_file2 = evec_file. \
                        replace('J'+str(j1), 'J'+str(j2)). \
                        replace('p'+str(p1), 'p'+str(p2))

                    path2 = os.path.join(evec_dir, evec_file2)
                    fcfr = -1

                    if os.path.isfile(path1) and os.path.isfile(path2):
                        fcfr = self.compute_single_FCF(
                            path1, path2, uch, lch, uv, lv,
                            shape, rgrid, igrid, istep
                        )

                        branch = self.map_parity_to_branch('R', p1, p2)
                        rbranch12 = np.array([
                            uv, lv, j1, j2, uch, lch, p1, p2, iso, fcfr, branch
                        ])

                        fcfs = np.hstack((fcfs, rbranch12))
                        counter += 1

                    ''' Pff/Pee branch '''

                    j2 = j1 + 1
                    p1 = int(evec_file.rsplit('_')[2][1:])
                    p2 = p1

                    # lower wavefunction
                    evec_file2 = evec_file. \
                        replace('J'+str(j1), 'J'+str(j2)). \
                        replace('p'+str(p1), 'p'+str(p2))

                    path2 = os.path.join(evec_dir, evec_file2)
                    fcfp = -1

                    if os.path.isfile(path1) and os.path.isfile(path2):
                        fcfp = self.compute_single_FCF(
                            path1, path2, uch, lch, uv, lv,
                            shape, rgrid, igrid, istep
                        )

                        branch = self.map_parity_to_branch('P', p1, p2)
                        pbranch12 = np.array([
                            uv, lv, j1, j2, uch, lch, p1, p2, iso, fcfp, branch
                        ])

                        fcfs = np.hstack((fcfs, pbranch12))
                        counter += 1

        return fcfs.reshape(counter, 11)

    def map_parity_to_branch(self, main_branch, p1, p2):

        qmap = {(1, 0): 'Qef', (0, 1): 'Qfe'}
        pmap = {(1, 1): 'Pee', (0, 0): 'Pff'}
        rmap = {(1, 1): 'Ree', (0, 0): 'Rff'}

        if main_branch == 'Q':
            return qmap[(p1, p2)]
        if main_branch == 'P':
            return pmap[(p1, p2)]
        if main_branch == 'R':
            return rmap[p1, p2]

    def compute_single_FCF(self, path1, path2, uch, lch, uv, lv,
                           shape, rgrid, igrid, istep):

        evectors1 = np.loadtxt(path1)
        evectors2 = np.loadtxt(path2)

        vibnums1 = evectors1[0, :].astype(np.int64)
        states1 = evectors1[1, :].astype(np.int64)

        vibnums2 = evectors2[0, :].astype(np.int64)
        states2 = evectors2[1, :].astype(np.int64)

        # find indices of upper and lower v's and states
        uvi, lvi = np.where(vibnums1 == uv)[0], np.where(vibnums2 == lv)[0]
        uchi, lchi = np.where(states1 == uch)[0], np.where(states2 == lch)[0]

        uevector = evectors1[2:, np.intersect1d(uvi, uchi)[0]]
        levector = evectors2[2:, np.intersect1d(lvi, lchi)[0]]

        # should be equal to the number of channels
        unch = int(uevector.shape[0] / shape)
        lnch = int(levector.shape[0] / shape)

        # start = time.time()

        # interpolate the eigenvector for the upper state
        upsi = self.interpolate_wavefunction(
            uevector, unch, rgrid, igrid, istep
        )
        upsi_norm = self.normilize_wavefunction(upsi, istep)

        # interpolate the eigenvector for the lower state
        lpsi = self.interpolate_wavefunction(
            levector, lnch, rgrid, igrid, istep
        )
        lpsi_norm = self.normilize_wavefunction(lpsi, istep)

        # print(f'wavefunc = {time.time()-start} s')

        overlap = np.dot(upsi_norm, lpsi_norm) * istep

        fcf = overlap * overlap  # * 1000.0

        return fcf

    def interpolate_wavefunction(self, evec, nch, rgrid, igrid, istep):
        # TODO: numpy implementation
        arr = []
        ngrid = rgrid.shape[0]
        ninter = igrid.shape[0]

        for row in range(0, ninter):
            index = 0
            sum_arg = 0.0
            for j in range(0, nch*ngrid):
                if j % ngrid == 0:
                    index += 1

                if abs(igrid[row] - rgrid[j-index*ngrid]) < 1e-15:
                    sinc = 1.0
                else:
                    arg = \
                        (math.pi / istep) * (igrid[row] - rgrid[j-index*ngrid])
                    sinc = math.sin(arg) / arg

                sum_arg += evec[j] * sinc

            arr.append(sum_arg)

        return np.array(arr)

        # # numpy implementation:
        # rotj = np.float64(os.path.splitext(efile)[0].split('_')[1])
        # evec = np.loadtxt(os.path.join(self.evec_dir, efile))
        # dif_shapes = abs(inter_grid.shape[0] - self.rgrid.shape[0])
        # rgrid = np.append(self.rgrid, self.rgrid[0:dif_shapes])
        # (np.pi / inter_step) * (inter_grid - rgrid)
        # arg = (1.0 / inter_step) * (inter_grid - rgrid)
        # sinc = np.sinc(arg) #np.sin(arg) / arg # np.sinc(x) -> sin(pi*x)/pi*x
        # sinc = np.tile(sinc[:,np.newaxis], (1, 30))
        # wavefunc = sinc * eig_vec #np.sum(eig_vec, axis=0)[:,np.newaxis].T
        # # normalized wavefunction
        # norm_wavefunc = wavefunc / wavefunc.sum(axis=0,keepdims=1)

    def normilize_wavefunction(self, psi, istep):

        # store the sign of psi
        psi_sign = np.sign(psi)

        # psi square
        psi_square = psi ** 2

        # the norm of psi square
        psi2_norm = psi_square / (psi_square.sum(axis=0, keepdims=1) * istep)

        psi_norm = psi_sign * np.sqrt(psi2_norm)

        return psi_norm

    def save_computed_FCFs(self, fname, fcfs):

        names = [
            "v'", "v", "J'", "J", "state'", "state",
            "p'", "p", "isotope", "FCF", "branch"
        ]

        fmt = [
            '%3.1s', '%6.1s', '%7.5s', '%7.5s', '%7.1s', '%7.1s',
            '%7.1s', '%7.1s', '%7.1s', '%13.10s', '%6s'
        ]

        header = \
            ''.join(['{:>3}'.format(k) for k in names[:1]]) + \
            ''.join(['{:>6}'.format(k) for k in names[1:2]]) + \
            ''.join(['{:>9}'.format(k) for k in names[2:3]]) + \
            ''.join(['{:>8}'.format(k) for k in names[3:8]]) + \
            ''.join(['{:>10}'.format(k) for k in names[8:]])

        np.savetxt(fname, fcfs, header=header, comments='#', fmt=fmt)

    # Cubic spline interpolation -> might not work correctly
    # for i in range(0, nch):
    #     interpolate the first eigenvector for the upper/lower state
    #     ynorm1 = self.interpolate_and_norm_evector(
    #       evector1, rgrid, shape, igrid, i)
    #     interpolate the second eigenvector for the lower/upper state
    #     ynorm2 = self.interpolate_and_norm_evector(
    #       evector2, rgrid, shape, igrid, i)

    # def interpolate_and_norm_evector(self, evector, rgrid, shape, igrid, i):
    #     cs = CubicSpline(
    #         rgrid,
    #        evector[i*shape:(i+1)*shape], bc_type='natural', extrapolate=True
    #     )
    #     y = cs(igrid)
    #     ynorm = np.zeros_like(y)
    #     if np.all((y != 0.0)):
    #         ynorm = y / y.sum(axis=0, keepdims=1)
    #     return ynorm
