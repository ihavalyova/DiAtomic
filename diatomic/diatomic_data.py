from re import compile as _compile, search as _search, findall as _findall
from io import open as _open
import numpy as np
from collections import OrderedDict
from .data.atomic_database import AtomicDatabase
from .utils import Utils as _utils


__all__ = ['Diatomic', 'Basis', 'Potential', 'Parameters']


class Diatomic:

    def __init__(self, molecule, masses=None, niso=None,
                 referencej=None, refE=None):

        self.molecule = molecule
        self.masses = None
        self.atomic_db = []
        self.reduced_masses = []
        self.atomic_symbols = []
        self.atomic_mass_nums = []
        self.referencej = referencej

        if molecule:
            for mname in self.molecule:
                rmass = self._calculate_reduced_mass(mname) * _utils.C_massau
                self.reduced_masses.append(rmass)
            ni = len(self.molecule)
        elif masses:
            self.masses = [m * _utils.C_massau for m in masses]
            ni = len(masses)
        else:
            raise ValueError(
                "Error: either 'molecule' or 'masses' should be provided!")

        self.niso = niso or list(range(1, ni+1))
        self.refE = refE / _utils.C_hartree if refE else None
        self.exp_file = None
        self.exp_data = None
        self.wavens_data = None
        self.uterms_file = None
        self.uterms_data = None
        self.lterms_file = None
        self.lterms_data = None

    def _calculate_reduced_mass(self, molecule, database=None):

        found = _search(r'^(\d+)(\D+)(\d+)(\D+)\b', molecule.strip())

        err_msg = f'Error: molecule {molecule} is not in the correct format'

        if not found:
            raise ValueError(err_msg)

        molecule_data = found.groups()

        if len(molecule_data) != 4:
            raise ValueError(err_msg)

        atoms = (molecule_data[1].strip(), molecule_data[3].strip())
        mass_numbers = (int(molecule_data[0]), int(molecule_data[2]))

        db_name = AtomicDatabase.database
        atomic_db = database or db_name

        atom_mass1 = self._match_atomic_mass(
            atomic_db, atoms[0], mass_numbers[0])
        atom_mass2 = self._match_atomic_mass(
            atomic_db, atoms[1], mass_numbers[1])

        rmass = atom_mass1 * atom_mass2 / (atom_mass1 + atom_mass2)

        self.atomic_db.append((atom_mass1, atom_mass2))
        self.atomic_symbols.append((atoms[0], atoms[1]))
        self.atomic_mass_nums.append((mass_numbers[0], mass_numbers[1]))

        return rmass

    def _match_atomic_mass(self, atomic_db, symbol, mnumber):

        pattern = _compile(
            fr'Atomic\s+Symbol\s+=\s+{symbol}[\r\n]'
            fr'Mass\s+Number\s+=\s+{mnumber}'
            fr'[\r\n]Relative\s+Atomic\s+Mass\s+=\s+\d+\.\d+')

        atom_data = _findall(pattern, atomic_db)

        if len(atom_data) != 1:
            raise ValueError(
                f'Error: nothing has been found for the symbol {symbol}.')

        return float(atom_data[0].split('\n')[-1].split('=')[-1].strip())

    def set_data_parameters(self, fname):

        self.params_by_labels = {}
        labels_pcount = {}
        xparams, yparams, fparams = [], [], []
        label = ''
        with open(fname) as fstream:
            for line in fstream:
                if line.startswith('@'):
                    if label != '':
                        self.params_by_labels[label] = np.column_stack(
                            (xparams, yparams, fparams))
                        xparams, yparams, fparams = [], [], []
                    label = line.strip(' @:\t\n')
                    labels_pcount[label] = 0
                elif line.lstrip()[0].isdigit():
                    labels_pcount[label] += 1
                    sline = line.split('#')[0]  # remove comments
                    sline = sline.strip().split()
                    if len(sline) == 3:
                        xparams.append(float(sline[0]))
                        yparams.append(float(sline[1]))
                        fparams.append(int(sline[2]))
                    elif len(sline) == 2:
                        xparams.append(0.0)  # dummy value
                        yparams.append(float(sline[0]))
                        fparams.append(int(sline[1]))
                    else:
                        raise ValueError(f'Wrong number of columns in {fname}')
                else:
                    continue

            self.params_by_labels[label] = np.column_stack((
                xparams, yparams, fparams))

        self.labels_inds = {}
        vstart = 0
        self.params = np.zeros(sum(labels_pcount.values()))
        for key, value in labels_pcount.items():
            vend = vstart + value
            self.labels_inds[key] = (vstart, vend)
            self.params[vstart:vend] = self.params_by_labels[key][:, 1]
            vstart += value

    def get_observed_energies(self):

        return self.exp_data

    def set_observed_energies(self, fname, ndata=None, markers=None,
                              apply_filters=False):

        self.exp_file = fname

        try:
            data = np.genfromtxt(self.exp_file, comments='#', autostrip=True)
            ndata = ndata or data.shape[0]
        except Exception as e:
            raise SystemExit(e)

        self.exp_data = self._apply_filter_to_observed_energies(
            ndata, data, markers, apply_filters)

        if self.exp_data.shape[0] == 0 or self.exp_data.shape[1] != 8:
            raise ValueError(
                f'Zero #levels or wrong format of exp energies in file {fname}'
                f'{self.exp_data.shape[0]} x {self.exp_data.shape[1]}\n')

    def _apply_filter_to_observed_energies(self, ndata, data, markers,
                                           apply_filters):

        # if data contains only a single line
        if data.ndim == 1:
            data = data[np.newaxis, :]

        # change the total number of levels
        data = data[:ndata, :]

        # filter by marker
        if markers:
            marker_mask = np.in1d(data[:, 6], markers)
            data = data[marker_mask]

        if apply_filters:
            # filter by symmetry
            parity_mask = np.in1d(data[:, 4], self.symmetry)
            data = data[parity_mask]

            # filter by J
            rot_mask = np.in1d(data[:, 2], self.jqnumbers)
            data = data[rot_mask]

            # change numbering
            data[:, 0] = np.arange(1.0, data.shape[0]+1)

        return data

    def set_observed_transitions(self, fname, markers=None,
                                 ithreshold=None):
        """Read the observed transitions data in specific format

        Args:
            fname (str): The name of the input file containing the data
            markers (list, optional): The list of markers saying which
                transitions to be ignored during the calculations.
                Defaults to None.

        Raises:
            IOError: If the file does not exist or cannot be loaded.
            ValueError: If the format of the file is wrong or
                the file does not contain any data.

        Note:
            If list of markers is provided, the last column of the data file
            is taken as marker column. The data are then filtered and the
            markers column is removed at the end.
        """

        try:
            self.wavens_data = np.loadtxt(fname)
            # self.euterms = np.unique(self.H.wavens_data[:, :4], axis=0)
            # self.elterms = np.unique(self.H.wavens_data[:, 4:8], axis=0)
        except IOError:
            raise IOError(f'Error: The file {fname} does not exist.')

        # if the data contains one line
        if self.wavens_data.ndim == 1:
            self.wavens_data = self.wavens_data[np.newaxis, :]

        rshape = self.wavens_data.shape[0]
        cshape = self.wavens_data.shape[1]

        ncols_no_mark = 12
        ncols_mark = 13

        if rshape == 0 or (cshape != ncols_no_mark and cshape != ncols_mark):
            raise ValueError(
                f'Cannot read the format of the observed transition. '
                f'Wrong shape of the data: {rshape} x {cshape} '
                f'or empty list.\n')

        # filter by marker
        if markers is not None:
            marker_mask = np.in1d(
                self.wavens_data[:, -1], np.fromiter(markers, dtype=np.int64))
            self.wavens_data = self.wavens_data[marker_mask]

            # remove the markers column
            self.wavens_data = self.wavens_data[:, :-1]

        # filter by intensity threshold
        if ithreshold is not None:
            max_intensity = np.amax(self.wavens_data[:, -2])
            ithresh_value = (ithreshold * max_intensity) / 100
            ithresh_inds = np.where(self.wavens_data[:, -2] > ithresh_value)[0]
            self.wavens_data = self.wavens_data[ithresh_inds, :]

    def set_uppper_term_energies(self, fname, markers=None,
                                 apply_filters=False):

        self.uterms_file = fname

        try:
            self.uterms_data = self._read_observed_energies(
                markers, apply_filters)
        except IOError:
            raise IOError(f'Error: The file {fname} does not exist.')

        if self.uterms_data.shape[0] == 0 or self.uterms_data.shape[1] != 8:
            raise ValueError(
                f'Wrong shape/format of observed energies in file {fname}'
                f'{self.uterms_data.shape[0]} x {self.uterms_data.shape[1]}\n'
                f'Check markers, Js, the format of the data and so on...')

    def set_lower_term_energies(self, fname, markers=None,
                                apply_filters=False):

        self.lterms_file = fname

        try:
            self.lterms_data = self._read_observed_energies(
                markers, apply_filters)
        except IOError:
            raise IOError(f'Error: The file {fname} does not exist.')

        if self.lterms_data.shape[0] == 0 or self.lterms_data.shape[1] != 8:
            raise ValueError(
                f'Wrong shape/format of observed energies in file {fname}'
                f'{self.lterms_data.shape[0]} x {self.lterms_data.shape[1]}\n'
                f'Check markers, Js, the format of the data and so on...')


class Basis:

    def __init__(self, nstate, Js=None, Je=None, Jvalues=None, symmetry=(0,),
                 _lambda=0, spin=0, sigma=0.0, inuc1=None, inuc2=None, F=None):

        self.nstate = nstate

        if Js is None and Je is None and Jvalues is None:
            raise ValueError('Error: Missing range/values for J quantum number')

        if Js is not None and Je is None:
            Je = Js + 1.0

        if Js is None and Je is not None:
            Js = Je - 1.0

        if Js < 0.0 or Je < 0.0:
            raise ValueError('Error: J quantum number should be positive')

        self.jqnumbers = []
        self.jqnumbers = np.arange(Js, Je, dtype=np.float64)

        if Jvalues is not None:
            if any(j < 0 for j in Jvalues):
                raise ValueError('Error: J quantum number should be positive')

            self.jqnumbers = np.unique(np.hstack((self.jqnumbers, Jvalues)))

        self.symmetry = (symmetry,) if isinstance(symmetry, int) else symmetry
        self._lambda = _lambda
        self.spin = spin
        self.sigma = sigma
        self.omega = abs(self._lambda + self.sigma)
        self.mult = 2 * self.spin + 1
        self.inuc1 = inuc1
        self.inuc2 = inuc2
        self.ftot = F

    def _format_basis_line(self, j, symm):

        maps = {0: 'f', 1: 'e'}
        out = (f'| \u039B={self._lambda}, S={self.spin}, \u03A3={self.sigma},'
               f' J={j}, \u03A9={self.omega}, symm={maps[symm]} ')

        if None not in [self.inuc1, self.inuc2, self.ftot]:
            out += f'I1={self.inuc1}, I2={self.inuc2}, F={self.ftot}'
        out += '>; '

        return out

    def print_basis_states(self):

        for j in self.jqnumbers:
            basis_str = self._format_basis_line(j, self.symmetry[0])
            try:
                basis_str += self._format_basis_line(j, self.symmetry[1])
            except IndexError:
                pass
            print(basis_str)


class Potential():

    def __init__(self, model, filep, rotc=0.0, shift_by=0.0, cfunc=None):

        self.model = model.lower()
        self.filep = filep
        self.rot_correction = rotc
        self.shift_by = shift_by
        self.cfunc = cfunc
        self.xunits = 1.0 / _utils.C_bohr
        self.yunits = 1.0 / _utils.C_hartree
        self.npnts = 0
        self.rpoints = None
        self.upoints = None
        self.fixed = None

        # regularization parameters
        self.pregular = 0
        self.plambda = 0
        self.id = id(self)

        if self.model == 'pointwise' or self.model == 'cspline':
            rpoints, upoints, fixed = self._read_pointwise_data(self.filep)
            self.xunits = np.full(rpoints.shape[0], self.xunits)
            self.yunits = np.full(upoints.shape[0], self.yunits)
            self.rpoints = rpoints
            self.upoints = upoints + self.shift_by
            self.fixed = fixed

        elif self.model == 'morse':
            self.upoints, self.fixed = self._read_morse_data(self.filep)

        elif self.model == 'emo':
            self.upoints, self.fixed = self._read_emo_data(self.filep)

        elif self.model == 'mlr':
            self.upoints, self.fixed = self._read_mlr_data(self.filep)

        elif self.model == 'custom':
            try:
                self.upoints, self.fixed = self._read_custom_data(self.filep)
            except TypeError as te:
                print(te)
        else:
            raise ValueError(f'Invalid potential model {self.model}')

        self.npnts = self.upoints.shape[0]

    def _read_pointwise_data(self, filep):

        points = np.loadtxt(filep)
        rpoints, upoints = points[:, 0], points[:, 1]

        fixed = np.zeros_like(upoints)
        if points.shape[1] == 3:
            fixed = points[:, 2]

        return rpoints, upoints, fixed

    def _read_morse_data(self, filep):

        # TODO: ignore comments and empty lines
        with open(filep) as fps:
            morse_data = fps.read().strip().split('\n')

        morse_data = dict(map(str.strip, s.split('=')) for s in morse_data)
        npt = len(morse_data)
        for md in morse_data.items():
            morse_data[md[0]] = md[1].split()

        def_fixed = False
        if len(list(morse_data.values())[0]) >= 2:
            def_fixed = True

        # TODO: check for key error

        mapp = {0: 'Te', 1: 'De', 2: 'a', 3: 're'}
        upoints, fixed = np.zeros(npt), np.zeros(npt, dtype=np.int64)

        upoints[0] = float(morse_data[mapp[0]][0]) / _utils.C_hartree
        upoints[1] = float(morse_data[mapp[1]][0]) / _utils.C_hartree
        upoints[2] = float(morse_data[mapp[2]][0]) * _utils.C_bohr
        upoints[3] = float(morse_data[mapp[3]][0]) / _utils.C_bohr

        if def_fixed:
            for i in range(0, 4):
                fixed[i] = int(morse_data[mapp[i]][1])

        return upoints, fixed

    def _read_emo_data(self, filep):

        # TODO: ignore comments and empty lines

        with open(filep) as fps:
            emo_data = fps.read().strip().split('\n')

        emo_data = dict(map(str.strip, s.split('=')) for s in emo_data)
        npt = len(emo_data)
        for md in emo_data.items():
            emo_data[md[0]] = md[1].split()

        def_fixed = False
        if len(list(emo_data.values())[0]) >= 2:
            def_fixed = True

        bparams = dict(
            filter(lambda x: x[0].lower().startswith('b'), emo_data.items()))

        # TODO: check for key error

        mapp = {0: 'Te', 1: 'De', 2: 'p', 3: 're'}

        upoints, fixed = np.zeros(npt), np.zeros(npt, dtype=np.int64)

        upoints[0] = float(emo_data[mapp[0]][0]) / _utils.C_hartree
        upoints[1] = float(emo_data[mapp[1]][0]) / _utils.C_hartree
        upoints[2] = float(emo_data[mapp[2]][0])
        upoints[3] = float(emo_data[mapp[3]][0]) / _utils.C_bohr

        # all beta coefficients have dimentions 1/distance
        bvalues = list(
            map(lambda x: float(x[0]) * _utils.C_bohr, bparams.values()))
        ni, nb = 4, len(bvalues)

        upoints[ni:ni+nb] = bvalues

        if def_fixed:
            bfixed = list(map(lambda x: int(x[1]), bparams.values()))

            for i in range(0, ni):
                fixed[i] = int(emo_data[mapp[i]][1])

            fixed[ni:ni+nb] = bfixed

        return upoints, fixed

    def _read_mlr_data(self, filep):
        """Get parameters for Morse/Long-Range function.

        Args:
            filep (string): the name of the potential file

        Returns:
            tuple of numpy arrays: the parameter values and free/fixed values

        Remarks:
            1. The coefficients B, C, D are defined from the first up to the
            highest reqired one. The ignored coefficients must be set to zero.
            2. The number of C and D coefficients should be the same
            3. The beta coefficients have dimentions 1/distance
        """
        # TODO: to ignore comments and empty lines

        with open(filep) as fps:
            mlr_data = fps.read().strip().split('\n')

        mlr_data = dict(map(str.strip, s.split('=')) for s in mlr_data)
        npt = len(mlr_data)
        for md in mlr_data.items():
            mlr_data[md[0]] = md[1].split()

        def_fixed = False
        if len(list(mlr_data.values())[0]) >= 2:
            def_fixed = True

        bparams = dict(
            filter(lambda x: x[0].lower().startswith('b'), mlr_data.items()))
        cparams = dict(
            filter(lambda x: x[0].lower().startswith('c'), mlr_data.items()))
        dparams = dict(
            filter(lambda x: x[0].lower().startswith('d'), mlr_data.items()))
        # remove De
        dparams = dict(
            filter(lambda x: x[0].lower() != 'de', dparams.items()))
        # remove binf
        bparams = dict(
            filter(lambda x: x[0].lower() != 'binf', bparams.items()))
        # TODO: check for key error

        mapp = {
            0: 'Te', 1: 'De', 2: 'p', 3: 'q',
            4: 'rref', 5: 're', 6: 'binf'
        }

        upoints, fixed = np.zeros(npt), np.zeros(npt, dtype=np.int64)
        upoints[0] = float(mlr_data[mapp[0]][0]) / _utils.C_hartree
        upoints[1] = float(mlr_data[mapp[1]][0]) / _utils.C_hartree
        upoints[2] = float(mlr_data[mapp[2]][0])
        upoints[3] = float(mlr_data[mapp[3]][0])
        upoints[4] = float(mlr_data[mapp[4]][0]) / _utils.C_bohr
        upoints[5] = float(mlr_data[mapp[5]][0]) / _utils.C_bohr
        upoints[6] = float(mlr_data[mapp[6]][0]) * _utils.C_bohr

        # TODO: check dimenstions of C and D parameters - they are not correct!
        bvalues = list(
            map(lambda x: float(x[0]) * _utils.C_bohr, bparams.values()))
        cvalues = list(
            map(lambda x: float(x[0]), cparams.values()))
        dvalues = list(
            map(lambda x: float(x[0]), dparams.values()))

        ni, nb, nc, nd = 7, len(bvalues), len(cvalues), len(dvalues)
        upoints[ni:ni+nb] = bvalues
        upoints[ni+nb:ni+nb+nc] = cvalues
        upoints[ni+nb+nc:ni+nb+nc+nd] = dvalues

        if def_fixed:
            bfixed = list(map(lambda x: int(x[1]), bparams.values()))
            cfixed = list(map(lambda x: int(x[1]), bparams.values()))
            dfixed = list(map(lambda x: int(x[1]), bparams.values()))

            for i in range(0, ni):
                fixed[i] = int(mlr_data[mapp[i]][1])

            fixed[ni:ni+nb] = bfixed
            fixed[ni+nb:ni+nb+nc] = cfixed
            fixed[ni+nb+nc:ni+nb+nc+nd] = dfixed

        # channel.npnts = ni + nb + nc + nd
        return upoints, fixed

    def _read_custom_data(self, filep):
        """Get parameters for a custom potential function.

        Args:
            filep (string): The name of the potential file

        Returns:
            tuple of two numpy arrays:
                the values of the parameters and the free/fixed values

        Remarks:
            1. The custom function should accept exactly 2 input parameters.
            The first one is an array containing the parameters and the
            second one is the grid points and should return one parameter-
            the values of the function on the grid points.
            3. The parameters should be in au units!
            4. Returned array has size equal to the number of grid points.
            5. A column for free/fixed value should alywas be defined
            6. The order in which the parameters are defined does not matter-
            they will be ordered by the number after the keyword 'param'!
        """
        # TODO: ignore comments and empty lines

        with open(filep, 'r') as fps:
            data = fps.read().strip().split('\n')

        data = OrderedDict(map(str.strip, s.split('=')) for s in data)
        data = OrderedDict(sorted(data.items()))

        for md in data.items():
            data[md[0]] = md[1].split()

        npt = len(data)

        upoints, fixed = np.zeros(npt), np.zeros(npt, dtype=np.int64)

        for i, v in enumerate(data.values()):
            upoints[i] = float(v[0])
            fixed[i] = int(v[1])
        
        return upoints, fixed

        # rpp, upoints, fixed = np.loadtxt(filep)
        # return upoints, fixed


class Channel:

    def __init__(self, id, basis, potential):

        self.basis = basis
        self.potential = potential
        self.upoints = self.potential.upoints
        self.fixed = self.potential.fixed
        self.cfunc = self.potential.cfunc

        # channels with the same poten file have the same id
        self.id = id
        self.refs = None

        # regularization parameters
        self.pregular = 0
        self.plambda = 0

    def write_pointwise_data(self, file_name, x, y, z):

        data = np.column_stack((x, y, z))
        fmt = ['%20.10f', '%25.10f', '%6d']
        np.savetxt(file_name, data, fmt=fmt)

    def write_morse_data(self, file_name, y, z):

        with open(file_name, 'w') as fs:
            fs.write(f'Te = {y[0]:>17.8e}{z[0]:>9}\n')
            fs.write(f'De = {y[1]:>17.8e}{z[1]:>9}\n')
            fs.write(f'a  = {y[2]:>17.8f}{z[2]:>9}\n')
            fs.write(f're = {y[3]:>17.8f}{z[3]:>9}\n')

    def write_emo_data(self, file_name, y, z):
        pass

    def write_mlr_data(self, file_name, y, z):
        pass

    def write_custom_data(self, file_name, y, z):

        with open(file_name, 'w') as f:
            for i in range(0, len(y)):
                f.write(f'param{i+1} = {y[i]:>17.8e}'
                        f'{z[i]:>10}\n')

    @classmethod
    def edit_channel_parameters(cls, ypar, channels):

        onpnts = 0

        for i, c in enumerate(channels):
            if i in c.id:
                npnts = c.upoints.shape[0]
                if c.model == 'pointwise' or c.model == 'cspline':
                    st, en = onpnts, onpnts + npnts
                    x = c.rpoints*_utils.C_bohr
                    y = ypar[st:en]*_utils.C_hartree
                    z = c.fixed
                    c.write_pointwise_data(c.filep, x, y, z)
                    onpnts = en
                elif c.model == 'morse':
                    st, en = onpnts, onpnts + npnts
                    Te, De, a, re = ypar[st:en]
                    y = [Te*_utils.C_hartree, De*_utils.C_hartree,
                         a, re*_utils.C_bohr]
                    c.write_morse_data(c.filep, y, c.fixed)
                    onpnts = en
                elif c.model == 'emo':
                    pass
                elif c.model == 'mlr':
                    pass
                elif c.model == 'custom':
                    st, en = onpnts, onpnts + npnts
                    params = ypar[st:en]
                    c.write_custom_data(c.filep, params, c.fixed)
                    onpnts = en


class Parameters:
    """A class setting and processing the fitting parameters
    """

    def __init__(self, params):

        # self.set_channels_refs(channels)
        # channel_pars = self.merge_channel_parameters(channels)

        # self.set_couplings_refs(couplings)
        # coupling_pars = self.merge_coupling_parameters(couplings)

        self.set_initial_parameters(channel_pars, coupling_pars)

        # total number of channels and coupling parameters
        self.tot_pars = self.tot_chnls_params + self.tot_cpls_params

    def set_initial_parameters(self, channel_pars, coupling_pars):

        self.ypar_init = np.concatenate((channel_pars[0], coupling_pars[0]))
        self.yfixed_init = np.concatenate((channel_pars[1], coupling_pars[1]))
        self.xunits = np.concatenate((channel_pars[2], coupling_pars[2]))
        self.yunits = np.concatenate((channel_pars[3], coupling_pars[3]))
        self.refs = np.concatenate((channel_pars[4], coupling_pars[4]))
        # self.refs = channel_pars[4]

        # this array will change during the fit
        self.ypar = np.copy(self.ypar_init)
        self.yfixed = np.copy(self.yfixed_init)

    def set_channels_refs(self, channels):

        ids = {}
        sind, eind = 0, 0
        for ch in channels:
            eind += ch.potential.npnts
            if ch.potential.id in ids:
                ch.refs = ids[ch.potential.id].refs
                ch.fixed = np.zeros(ch.potential.npnts)
            else:
                ch.refs = np.arange(sind, eind)
                ids[ch.potential.id] = ch
            sind = eind

    def set_couplings_refs(self, couplings):

        sind, eind = self.tot_chnls_params, self.tot_chnls_params

        for cp in couplings:
            cp.refs = np.zeros(len(cp.interact)*cp.npnts, dtype=int)
            cp.refs[0:cp.npnts] = np.arange(sind, sind+cp.npnts)
            cp.yfixed = np.zeros(len(cp.interact)*cp.npnts, dtype=int)
            cp.yfixed[0:cp.fixed.shape[0]] = cp.fixed
            eind += cp.npnts

            for i in range(1, len(cp.interact)):
                cp.refs[i*cp.npnts:(i+1)*cp.npnts] = np.arange(sind, eind)
                cp.yfixed[i*cp.npnts:(i+1)*cp.npnts] = np.zeros(cp.npnts, dtype=int)

            sind += len(cp.interact)*cp.npnts
            eind = sind

    def merge_channel_parameters(self, channels):

        ppar = np.concatenate([c.upoints for c in channels])
        pfixed = np.concatenate([c.fixed for c in channels])

        pxunits = np.concatenate([c.potential.xunits for c in channels])
        pyunits = np.concatenate([c.potential.yunits for c in channels])

        refs = np.concatenate([c.refs for c in channels])
        # pregular = np.concatenate([c.pregular for c in channels])
        # plambda = np.concatenate([c.plambda for c in channels])

        self.tot_chnls_params = ppar.shape[0]

        return ppar, pfixed, pxunits, pyunits, refs

    def merge_coupling_parameters(self, couplings):

        self.tot_cpls_params = 0

        try:
            cpar = np.concatenate(
                [c.yc for c in couplings for i in c.interact])
            cfixed = np.concatenate(
                [c.yfixed for c in couplings])
            cxunits = np.concatenate(
                [c.xunits for c in couplings for i in c.interact])
            cyunits = np.concatenate(
                [c.yunits for c in couplings for i in c.interact])
            crefs = np.concatenate(
                [c.refs for c in couplings])
            # self.cregular = np.concatenate([c.cregular for c in couplings])
            # self.clambda = np.concatenate([c.clambda for c in couplings])
            self.tot_cpls_params = cpar.shape[0]
        except ValueError:
            cpar, cfixed = np.array([]), np.array([], dtype=np.int64)
            cxunits, cyunits = np.array([]), np.array([])
            crefs = np.array([], dtype=np.int64)

        return cpar, cfixed, cxunits, cyunits, crefs

    def set_dmf_parameters(self, dmfs_init):

        self.tot_dmfs_params = self.tot_pars

        for (n, k) in dmfs_init:
            ypar_shape_old = self.ypar.shape[0]
            dmf_shape = dmfs_init[(n, k)].shape[0]

            self.ypar_init = np.concatenate(
                (self.ypar_init, dmfs_init[(n, k)][:, 1]))

            self.yfixed_init = np.concatenate(
                (self.yfixed_init, dmfs_init[(n, k)][:, 2]))

            self.ypar = np.concatenate((self.ypar, dmfs_init[(n, k)][:, 1]))

            ypar_shape_new = self.ypar.shape[0]

            self.yfixed = np.concatenate((self.ypar, dmfs_init[(n, k)][:, 2]))
            self.xunits = np.concatenate((self.xunits, np.full(dmf_shape, 1)))
            self.yunits = np.concatenate((self.yunits, np.full(dmf_shape, 1)))

            self.refs = np.concatenate(
                (self.refs, np.arange(ypar_shape_old, ypar_shape_new)))

            self.tot_dmfs_params += dmf_shape


class Coupling:

    def __init__(self, interact, coupling, model, label, multiplier=1,
                 custom_func=None):

        self.interact = interact
        self.coupling = coupling
        self.model = model.lower()
        self.label = label
        self.multiplier = multiplier
        self.cfunc = custom_func
        self.xunits = 1.0 / _utils.C_bohr
        self.yunits = 1
        self.npnts = 0
        self.fixed = 0
        self.xc = None
        self.yc = None
        self.refs = None

        # regularization parameters
        self.cregular = 0
        self.clambda = 0

        # parse input parameters
        parstr = list(map(str.split, Coupling.cpl_data[self.label]))
        params = np.array([list(map(float, i)) for i in parstr])

        # if string convert to tuples
        if isinstance(self.coupling, str):
            self.coupling = (self.coupling,)
            self.interact = (self.interact,)
            self.multiplier = (self.multiplier,)

        # convert to lowercase
        self.coupling = tuple(map(str.lower, self.coupling))

        # set units
        if 'spin-orbit' in self.coupling or 'dbobc' in self.coupling:
            self.yunits = 1.0 / _utils.C_hartree

        if all(item.startswith('lambdad') for item in self.coupling):
            self.yunits = _utils.C_hartree

        # set parameters by checking the model
        if self.model == 'pointwise' or self.model == 'cspline':
            self.xunits = np.full(params[:, 0].shape[0], self.xunits)
            self.yunits = np.full(params[:, 1].shape[0], self.yunits)
            self.xc = params[:, 0]
            self.yc = params[:, 1]

            if params.shape[1] >= 3:
                self.fixed = params[:, 2]
            if params.shape[1] >= 4:
                self.cregular = params[:, 3]
                self.clambda = params[:, 4]

        if self.model == 'custom':
            self.yc = params[:, 0]
            if params.shape[1] >= 2:
                self.fixed = params[:, 1]
            if params.shape[1] >= 3:
                self.cregular = params[:, 2]
                self.clambda = params[:, 3]

        self.npnts = params.shape[0]

    # @classmethod
    # def set_couplings_file(cls, fname):

    #     cls.cpl_file = fname

    #     try:
    #         with open(cls.cpl_file, 'r') as inps:
    #             try:
    #                 cls.cpl_data = yaml.safe_load(inps)
    #             except yaml.YAMLError as ye:
    #                 # TODO: fails with YAML object has no attribute YAMLError
    #                 raise SystemExit(ye)
    #     except IOError as e:
    #         raise SystemExit(e)

    @classmethod
    def edit_coupling_parameters(cls, ypar, couplings):

        cpl_data = Diatomic._read_couplings_data()
        # cpar = tot_npts
        cpar = 1

        for cp in couplings:
            new_item = []
            if cp.model == 'pointwise' or cp.model == 'cspline':
                for item in cpl_data[cp.label]:
                    sitem = item.split()
                    if len(sitem) >= 3:
                        new_item.append(
                            f'{float(sitem[0]):10.12f}'
                            f'{float(ypar[cpar]/cp.yunits):24.12f}'
                            f'{int(sitem[2]):7d}'
                            f'{float(sitem[3]):14.2e}'
                            f'{float(sitem[4]):14.1e}'.lstrip())
                    else:
                        new_item.append(
                            f'{float(sitem[0]):10.12f}'
                            f'{float(ypar[cpar]/cp.yunits):24.12f}'
                            f'{int(sitem[2]):7d}')
                    cpar += 1
            else:
                for item in cpl_data[cp.label]:
                    sitem = item.split()
                    if len(sitem) >= 2:
                        new_item.append(
                            f'{float(ypar[cpar]):18.12f}'
                            f'{int(sitem[1]):3d}'
                            f'{float(sitem[2]):10.1f}'.lstrip())
                    else:
                        new_item.append(
                            f'{float(ypar[cpar]):18.12f}'
                            f'{int(sitem[1]):3d}'.lstrip())
                    cpar += 1
            cpl_data[cp.label] = new_item

        # with _open(Diatomic.cpl_file, 'w', encoding='utf8') as stream:
        #     yaml.dump(cpl_data, stream=stream)
