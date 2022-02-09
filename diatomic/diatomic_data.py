from re import compile as _compile, search as _search, findall as _findall
from io import open as _open
import numpy as np
from collections import OrderedDict
from data.atomic_database import AtomicDatabase
from utils import C_hartree, C_bohr, C_massau

# to check if C loader is present on current machine
# /path/to/python -c "import ruamel.yaml;
# print(ruamel.yaml.__with_libyaml__)"

from ruamel.yaml import YAML
yaml = YAML(typ='rt', pure=False)

__all__ = ['DiatomicData']


class DiatomicData:

    def __init__(self, molecule, masses=None, niso=None, jrange=None,
                 jvalues=None, symmetry=(0, 1), refJ=None, refE=None):

        self.atomic_db = []
        self.reduced_masses = []
        self.masses = None
        self.atomic_symbols = []
        self.atomic_mass_nums = []
        self.molecule = molecule

        if molecule:
            for molecule_name in self.molecule:
                rmass = self._calculate_reduced_mass(molecule_name) * C_massau
                self.reduced_masses.append(rmass)
            ni = len(self.molecule)

        elif masses:
            self.masses = [m * C_massau for m in masses]
            ni = len(masses)

        else:
            raise ValueError(
                "One of the options 'molecule' or 'masses' is required!")

        self.niso = niso or list(range(1, ni+1))

        self.jqnumbers = []

        if jrange:
            self.jqnumbers = np.arange(jrange[0], jrange[1], dtype=np.float64)

        if jvalues:
            self.jqnumbers = np.unique(np.hstack((self.jqnumbers, jvalues)))

        self.symmetry = (symmetry,) if isinstance(symmetry, int) else symmetry

        self.refE = refE / C_hartree if refE else None
        self.refj = refJ or None
        if self.refj:
            self._arange_jqnumbers()

        self.exp_file = None
        self.exp_data = None
        self.wavens_data = None
        self.uterms_file = None
        self.uterms_data = None
        self.lterms_file = None
        self.lterms_data = None

    def _arange_jqnumbers(self):

        self.jqnumbers = np.delete(
            self.jqnumbers, np.argwhere(self.jqnumbers == self.refj))

        self.jqnumbers = np.insert(self.jqnumbers, 0, self.refj)

    def _calculate_reduced_mass(self, molecule, database=None):

        found = _search(r'^(\d+)(\D+)(\d+)(\D+)\b', molecule.strip())

        err_msg = f'Error: Molecule {molecule} not in the correct format'

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
                f'Error: Incorrect matching or nothing found in '
                f'the atomic database for symbol {symbol}.')

        return float(atom_data[0].split('\n')[-1].split('=')[-1].strip())

    def get_observed_energies(self):
        return self.exp_data

    def set_observed_energies(self, file_name, markers=None,
                              apply_filters=False):

        self.exp_file = file_name

        try:
            self.exp_data = self._read_observed_energies(
                markers, apply_filters)
        except Exception as e:
            raise SystemExit(e)

        if self.exp_data.shape[0] == 0 or self.exp_data.shape[1] != 8:
            raise ValueError(
                f'Wrong shape/format of observed energies in file {file_name}'
                f'{self.exp_data.shape[0]} x {self.exp_data.shape[1]}\n')

    def _read_observed_energies(self, markers, apply_filters):

        with open(self.exp_file) as efile:
            ndata = int(efile.readline())

        data = np.genfromtxt(self.exp_file, skip_header=1, comments='#',
                             max_rows=int(ndata), autostrip=True)

        # if the data contains only a single line
        if data.ndim == 1:
            data = data[np.newaxis, :]

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

    def set_observed_wavenumbers(self, file_name, markers=None):

        try:
            self._read_observed_wavenumbers(file_name, markers)
            # self.euterms = np.unique(self.H.wavens_data[:, :4], axis=0)
            # self.elterms = np.unique(self.H.wavens_data[:, 4:8], axis=0)
        except Exception as err:
            raise SystemExit(err)

        ncols = 12

        if self.wavens_data.shape[0] == 0 or \
           self.wavens_data.shape[1] != ncols:
            raise ValueError(
                f'Cannot read the format of the observed wavenumbers. '
                f'Wrong shape of the data: '
                f'{self.wavens_data.shape[0]} x {self.wavens_data.shape[1]}\n')

    def _read_observed_wavenumbers(self, file_name, markers):

        self.wavens_data = np.loadtxt(file_name)

        # if the data contains one line
        if self.wavens_data.ndim == 1:
            self.wavens_data = self.wavens_data[np.newaxis, :]

        # filter by marker
        # if markers is not None:
        #     marker_mask = np.in1d(
        #         self.wavens_data[:, 6], np.fromiter(markers, dtype=np.int64))
        #     self.wavens_data = self.wavens_data[marker_mask]

    def set_uppper_term_energies(self, file_name, markers=None,
                                 apply_filters=False):

        self.uterms_file = file_name

        try:
            self.uterms_data = self._read_observed_energies(
                markers, apply_filters)
        except Exception as e:
            raise SystemExit(e)

        if self.uterms_data.shape[0] == 0 or self.uterms_data.shape[1] != 8:
            raise ValueError(
                f'Wrong shape/format of observed energies in file {file_name}'
                f'{self.uterms_data.shape[0]} x {self.uterms_data.shape[1]}\n'
                f'Check markers, Js, the format of the data and so forth...')

    def set_lower_term_energies(self, file_name, markers=None,
                                apply_filters=False):

        self.lterms_file = file_name

        try:
            self.lterms_data = self._read_observed_energies(
                markers, apply_filters)
        except Exception as e:
            raise SystemExit(e)

        if self.lterms_data.shape[0] == 0 or self.lterms_data.shape[1] != 8:
            raise ValueError(
                f'Wrong shape/format of observed energies in file {file_name}'
                f'{self.lterms_data.shape[0]} x {self.lterms_data.shape[1]}\n'
                f'Check markers, Js, the format of the data and so forth...')

    def get_channel_parameters(self, channels):

        self.unq_channels = OrderedDict()
        count_pnts = 0

        # keep here the indices of the unique channels
        self.unq_chind = []

        # TODO: check if counting is correct
        for i, ch in enumerate(channels):
            if ch.filep not in self.unq_channels:
                self.unq_channels[ch.filep] = (count_pnts, count_pnts+ch.npnts)
                count_pnts += ch.npnts
                self.unq_chind.append(i)

            ch.start_index = self.unq_channels[ch.filep][0]
            ch.end_index = self.unq_channels[ch.filep][1]

        ppar = np.concatenate(
            [c.upoints for i, c in enumerate(channels) if i in self.unq_chind])
        pfixed = np.concatenate(
            [c.fixed for i, c in enumerate(channels) if i in self.unq_chind])
        pxunits = np.concatenate(
            [c.xunits for i, c in enumerate(channels) if i in self.unq_chind])
        pyunits = np.concatenate(
            [c.yunits for i, c in enumerate(channels) if i in self.unq_chind])

        # pregular = np.concatenate([c.pregular for c in channels])
        # plambda = np.concatenate([c.plambda for c in channels])

        # the total number of potential parameters
        self.tot_npts = ppar.shape[0]

        return ppar, pfixed, pxunits, pyunits

    def get_coupling_parameters(self, couplings):

        # TODO: check if counting is correct
        count_pnts = self.tot_npts
        for i, cp in enumerate(couplings):
            cp.start_index = count_pnts
            count_pnts += cp.npnts
            cp.end_index = count_pnts

        try:
            cpar = np.concatenate([c.yc for c in couplings])
            cfixed = np.concatenate([c.fixed for c in couplings])
            cxunits = np.concatenate([c.xunits for c in couplings])
            cyunits = np.concatenate([c.yunits for c in couplings])
            # self.cregular = np.concatenate([c.cregular for c in couplings])
            # self.clambda = np.concatenate([c.clambda for c in couplings])
        except ValueError:
            cpar, cfixed = np.array([]), np.array([])
            cxunits, cyunits = np.array([]), np.array([])

        return cpar, cfixed, cxunits, cyunits

    @classmethod
    def set_couplings_data(cls, file_name):

        cls.cpl_file = file_name
        cls.cpl_data = cls.read_couplings_data()

    @classmethod
    def read_couplings_data(cls):

        try:
            with open(cls.cpl_file, 'r') as inps:
                try:
                    return yaml.load(inps)
                except yaml.YAMLError as exc:
                    # TODO: fails with YAML object has no attribute YAMLError
                    raise SystemExit(exc)
        except IOError as e:
            raise SystemExit(e)

    def edit_channel_parameters(self, ypar, channels):

        onpnts = 0

        for i, c in enumerate(channels):
            if i in self.unq_chind:
                npnts = c.upoints.shape[0]
                if c.model == 'pointwise' or c.model == 'cspline':
                    st, en = onpnts, onpnts + npnts
                    x, y, z = c.rpoints*C_bohr, ypar[st:en]*C_hartree, c.fixed
                    c.write_pointwise_data(c.filep, x, y, z)
                    onpnts = en
                elif c.model == 'morse':
                    st, en = onpnts, onpnts + npnts
                    Te, De, a, re = ypar[st:en]
                    y = [Te*C_hartree, De*C_hartree, a, re*C_bohr]
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

    def edit_coupling_parameters(self, ypar, couplings):

        cpl_data = DiatomicData._read_couplings_data()
        cpar = self.tot_npts

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

        with _open(DiatomicData.cpl_file, 'w', encoding='utf8') as stream:
            yaml.dump(cpl_data, stream=stream)


class Channel:

    def __init__(self, **kwargs):

        self.filep = kwargs['filep']
        self.model = kwargs['model'].lower()
        self.nlambda = kwargs['nlambda']
        self.nsigma = kwargs['sigma']
        self.omega = abs(self.nlambda + self.nsigma)
        self.spin = kwargs['spin']
        self.mult = 2*self.spin + 1
        self.rot_correction = kwargs.get('rotc') or 0.0
        self.custom_function = kwargs.get('custom_function') or None
        self.start_index = None  # start index in ypar
        self.end_index = None  # end index in ypar
        self.xunits = 1.0 / C_bohr
        self.yunits = 1.0 / C_hartree
        self.npnts = 0
        self.rpoints = None
        self.upoints = None
        self.fixed = 0
        self.is_unique = False

        # channels with the same potential file will have the same id
        self.id = 0

        # regularization parameters
        self.pregular = 0
        self.plambda = 0

        if self.model == 'pointwise' or self.model == 'cspline':
            rpoints, upoints, fixed = self._read_pointwise_data(self.filep)
            self.xunits = np.full(rpoints.shape[0], self.xunits)
            self.yunits = np.full(upoints.shape[0], self.yunits)
            self.rpoints = rpoints * self.xunits
            self.upoints = upoints * self.yunits
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
                self.cfunc = self.custom_function
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

        upoints[0] = float(morse_data[mapp[0]][0]) / C_hartree
        upoints[1] = float(morse_data[mapp[1]][0]) / C_hartree
        upoints[2] = float(morse_data[mapp[2]][0]) * C_bohr
        upoints[3] = float(morse_data[mapp[3]][0]) / C_bohr

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

        upoints[0] = float(emo_data[mapp[0]][0]) / C_hartree
        upoints[1] = float(emo_data[mapp[1]][0]) / C_hartree
        upoints[2] = float(emo_data[mapp[2]][0])
        upoints[3] = float(emo_data[mapp[3]][0]) / C_bohr

        # all beta coefficients have dimentions 1/distance
        bvalues = list(map(lambda x: float(x[0]) * C_bohr, bparams.values()))
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
        upoints[0] = float(mlr_data[mapp[0]][0]) / C_hartree
        upoints[1] = float(mlr_data[mapp[1]][0]) / C_hartree
        upoints[2] = float(mlr_data[mapp[2]][0])
        upoints[3] = float(mlr_data[mapp[3]][0])
        upoints[4] = float(mlr_data[mapp[4]][0]) / C_bohr
        upoints[5] = float(mlr_data[mapp[5]][0]) / C_bohr
        upoints[6] = float(mlr_data[mapp[6]][0]) * C_bohr

        # TODO: check dimenstions of C and D parameters - they are not correct!
        bvalues = list(
            map(lambda x: float(x[0]) * C_bohr, bparams.values()))
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


class Coupling:

    def __init__(self, **kwargs):

        self.interact = kwargs['interact']
        self.coupling = kwargs['coupling']
        self.model = kwargs['model'].lower()
        self.label = kwargs['label']
        self.multiplier = kwargs.get('multiplier') or 1.0
        self.custom_function = kwargs.get('custom_function') or None
        self.xc = None
        self.yc = None
        self.npnts = 0
        self.fixed = 0
        self.xunits = 1.0 / C_bohr
        self.yunits = 1
        self.start_index = None  # start index in ypar
        self.end_index = None  # end index in ypar

        # regularization parameters
        self.cregular = 0
        self.clambda = 0

        params_str = list(map(str.split, DiatomicData.cpl_data[self.label]))
        params = np.array([list(map(float, i)) for i in params_str])

        # if string convert to tuples
        if isinstance(self.coupling, str):
            self.coupling = (self.coupling,)
            self.interact = (self.interact,)
            self.multiplier = (self.multiplier,)

        # convert to lowercase
        self.coupling = tuple(map(str.lower, self.coupling))

        # set units
        if 'spin-orbit' in self.coupling or 'dbobc' in self.coupling:
            self.yunits = 1.0 / C_hartree

        if all(item.startswith('lambdad') for item in self.coupling):
            self.yunits = C_hartree

        if self.model == 'pointwise' or self.model == 'cspline':
            self.xc = params[:, 0] * self.xunits
            self.yc = params[:, 1] * self.yunits

            if params.shape[1] >= 3:
                self.fixed = params[:, 2]
            if params.shape[1] >= 4:
                self.cregular = params[:, 3]
                self.clambda = params[:, 4]

        if self.model == 'constant' or self.model == 'custom':
            self.yc = params[:, 0]
            if params.shape[1] >= 2:
                self.fixed = params[:, 1]
            if params.shape[1] >= 3:
                self.cregular = params[:, 2]
                self.clambda = params[:, 3]

        self.cfunc = self.custom_function
        self.npnts = params.shape[0]
