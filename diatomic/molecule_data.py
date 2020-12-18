from re import compile as _compile, search as _search, findall as _findall
from os.path import splitext as _splitext
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

__all__ = ['MoleculeData', 'Channel', 'Coupling']


class MoleculeData:

    def __init__(self):

        self.imasses = None
        self.molecule_names = None
        self.pars = (0, 1)
        self.niso = None
        self.rot_range = None
        self.rot_values = None
        self.refj = None
        self.refE = None
        self.exp_file = None
        self.exp_data = None
        self.atomic_db = []
        self.reduced_masses = []
        self.atomic_symbols = []
        self.atomic_mass_nums = []
        self.jqnumbers = np.array([], dtype=np.float64)

    @property
    def molecule(self):
        return self.molecule_names

    @molecule.setter
    def molecule(self, molecule_names):
        self.reduced_masses = []
        self.molecule_names = molecule_names

        for molecule_name in self.molecule_names:
            rmass = self._calculate_reduced_mass(molecule_name) * C_massau
            self.reduced_masses.append(rmass)

    @property
    def masses(self):
        return self.imasses

    @masses.setter
    def masses(self, imasses):
        self.imasses = [m * C_massau for m in imasses]

    @property
    def nisotopes(self):
        return self.niso

    @nisotopes.setter
    def nisotopes(self, niso):
        self.niso = niso

    @property
    def jrange(self):
        return self.jqnumbers

    @jrange.setter
    def jrange(self, rot_range):
        self.rot_range = rot_range
        self.jqnumbers = np.hstack((
            self.jqnumbers,
            np.arange(rot_range[0], rot_range[1]+1.0, dtype=np.float64)
        ))
        self.jqnumbers = np.unique(self.jqnumbers)

    @property
    def jvalues(self):
        return self.jqnumbers

    @jvalues.setter
    def jvalues(self, rot_values):
        self.rot_values = rot_values
        self.jqnumbers = np.unique(
            np.hstack((self.jqnumbers, rot_values))
        )

    @property
    def referencej(self):
        return self.refj

    @referencej.setter
    def referencej(self, refj):
        self.refj = refj
        self._arange_jqnumbers()

    @property
    def parities(self):
        return self.pars

    @parities.setter
    def parities(self, pars):
        if isinstance(pars, int):
            pars = pars,

        self.pars = pars

    @property
    def referenceE(self):

        return self.refE

    @referenceE.setter
    def referenceE(self, refE):

        self.refE = refE / C_hartree

    def _arange_jqnumbers(self):

        self.jqnumbers = np.delete(
            self.jqnumbers, np.argwhere(self.jqnumbers == self.refj)
        )
        self.jqnumbers = np.insert(self.jqnumbers, 0, self.refj)

    def _calculate_reduced_mass(self, molecule, database=None):

        found = _search(r'^(\d+)(\D+)(\d+)(\D+)\b', molecule.strip())

        err_msg = f'Error: Molecule {molecule} not in the correct format'

        if found is None:
            raise SystemExit(err_msg)

        molecule_data = found.groups()

        if len(molecule_data) != 4:
            raise SystemExit(err_msg)

        atoms = (molecule_data[1].strip(), molecule_data[3].strip())
        mass_numbers = (int(molecule_data[0]), int(molecule_data[2]))

        db_name = AtomicDatabase.database
        atomic_db = database or db_name

        atom_mass1 = self._match_atomic_mass(
            atomic_db, atoms[0], mass_numbers[0]
        )
        atom_mass2 = self._match_atomic_mass(
            atomic_db, atoms[1], mass_numbers[1]
        )

        rmass = atom_mass1 * atom_mass2 / (atom_mass1 + atom_mass2)

        self.atomic_db.append((atom_mass1, atom_mass2))
        self.atomic_symbols.append((atoms[0], atoms[1]))
        self.atomic_mass_nums.append((mass_numbers[0], mass_numbers[1]))

        return rmass

    def _match_atomic_mass(self, atomic_db, symbol, mnumber):

        pattern = _compile(
            fr'Atomic\s+Symbol\s+=\s+{symbol}[\r\n]'
            fr'Mass\s+Number\s+=\s+{mnumber}'
            fr'[\r\n]Relative\s+Atomic\s+Mass\s+=\s+\d+\.\d+'
        )

        atom_data = _findall(pattern, atomic_db)

        if len(atom_data) != 1:
            raise SystemExit(
                f'Error: Incorrect matching or nothing found in '
                f'the atomic database for symbol {symbol}.'
            )

        return float(atom_data[0].split('\n')[-1].split('=')[-1].strip())

    def get_exp_data(self):

        return self.exp_data

    def set_exp_data(self, exp_file, markers=None, average=False):

        self.exp_file = exp_file

        try:
            self._read_experimental_data(markers, average)
        except Exception as e:
            # _error('Failed to read the experimental data: '+ str(e))
            raise SystemExit(e)

        if self.exp_data.shape[0] == 0 or self.exp_data.shape[1] != 8:
            raise SystemExit(
                f'Wrong shape of experimental data: '
                f'{self.exp_data.shape[0]} x {self.exp_data.shape[1]}\n'
                f'Check markers, Js, the format of the data and so forth...'
            )

    def _read_experimental_data(self, markers, average):

        ndata = np.genfromtxt(
            self.exp_file, max_rows=1, comments='#', autostrip=True
        )

        self.exp_data = np.genfromtxt(
            self.exp_file, skip_header=1, max_rows=int(ndata),
            comments='#', autostrip=True
        )

        # if the data contains one line
        if self.exp_data.ndim == 1:
            self.exp_data = self.exp_data[np.newaxis, :]

        # filter by marker
        if markers is not None:
            marker_mask = np.in1d(
                self.exp_data[:, 6], np.fromiter(markers, dtype=np.int64)
            )
            self.exp_data = self.exp_data[marker_mask]

        # filter by parity
        # if len(self.pars) == 1:
        self.pars = np.intersect1d(
            self.pars, np.unique(self.exp_data[:, 4])
        )

        parity_mask = np.in1d(
            self.exp_data[:, 4], np.fromiter(self.pars, dtype=np.int64)
        )
        self.exp_data = self.exp_data[parity_mask]

        self.jqnumbers = np.intersect1d(
            self.jqnumbers, np.unique(self.exp_data[:, 2])
        )

        if self.refj is not None:
            self._arange_jqnumbers()

        # filter by J
        rot_mask = np.in1d(
            self.exp_data[:, 2], np.fromiter(self.jqnumbers, dtype=np.float64)
        )
        self.exp_data = self.exp_data[rot_mask]

        self.exp_data[:, 0] = np.arange(1.0, self.exp_data.shape[0]+1)

        if average:
            self._average_experimental_data()

    def _average_experimental_data(self):

        # TODO: will not work when markers are not provided
        # TODO: to account for different isotopes
        # TODO: add weighted average of the exp uncer.

        # sort experimental data by v, J, parity and state
        self.exp_data = self.exp_data[np.lexsort((
            self.exp_data[:, -1], self.exp_data[:, 4],
            self.exp_data[:, 2], self.exp_data[:, 1]
        ))]

        exp_extract = self.exp_data[:, [1, 2, 4, -1]]
        unique_data, inds, counts = np.unique(
            exp_extract, axis=0, return_index=True, return_counts=True
        )

        avg_data = self.exp_data[inds]

        # create a generator yielding the average values
        # rep_gen = (self.exp_data[~np.any(exp_extract-row, axis=1)] \
        #   for row in unique_data[counts>1])

        for row in (unique_data[counts > 1]):
            item = self.exp_data[~np.any(exp_extract-row, axis=1)]

            avg_item = np.concatenate((
                    item[0, [0, 1, 2]],
                    np.array([np.average(item[:, 3])]),
                    item[0, [4]],
                    np.array([np.average(item[:, 5])]),
                    item[0, [6, 7]]
                ),
                axis=0
            )[np.newaxis, :]

            mask = np.column_stack((
                avg_data[:, [1, 2, 4, 7]] == avg_item[:, [1, 2, 4, 7]],
                np.full((avg_data.shape[0], 4), True)
            ))

            mask = np.all(mask, axis=1)
            avg_data[mask] = avg_item

        avg_data[:, 0] = np.arange(1.0, avg_data.shape[0]+1)

        self._save_averaged_experimental_data(avg_data)

    def _save_averaged_experimental_data(self, avg_data):

        header = str(avg_data.shape[0]) + '\n' + '# markers: ' + \
            np.array2string(
                np.unique(self.exp_data[:, 6]),
                formatter={'float_kind': lambda x: "%d" % x}
            )[1:-1]

        fmt = 2*['%5d'] + ['%7.1f', '%15.6f', '%7d', '%11.5f'] + 2*['%6d']

        fname, fext = _splitext(self.exp_file)
        avg_file = fname + '_averaged' + fext

        np.savetxt(avg_file, avg_data, header=header, comments='', fmt=fmt)


class Channel:

    def __init__(self, **kwargs):

        self.filep = kwargs['filep']
        self.model = kwargs['model'].lower()
        self.nlambda = kwargs['nlambda']
        self.nsigma = kwargs['sigma']
        self.omega = abs(self.nlambda + self.nsigma)
        self.mult = kwargs['multiplicity']
        self.rot_correction = kwargs.get('rotc') or 0.0
        self.custom_function = kwargs.get('custom_function')

    @classmethod
    def get_channel_objects(cls):
        return cls.channels

    @classmethod
    def _define_channel_models(cls):

        return {
            1: 'pointwise',
            2: 'cspline',
            3: 'morse',
            4: 'emo',
            5: 'mlr',
            6: 'custom',
            7: 'cpe'
        }

    @classmethod
    def set_channel_parameters(cls, channels):

        cls.models = cls._define_channel_models()
        cls.channels = channels
        cls.unique_channels = OrderedDict()
        cls.totp, cls.totp_unique = 0, 0

        for ci, ch in enumerate(channels):

            if ch.model == cls.models[1] or ch.model == cls.models[2]:
                ch.yunits = 1.0 / C_hartree
                ch.R, Uy, ch.fixedU = cls._get_pointwise_data(ch.filep)
                ch.U = Uy * ch.yunits
                ch.npnts = ch.U.shape[0]

            elif ch.model == cls.models[3]:
                ch.U, ch.fixedU = cls._get_morse_data(ch.filep)
                ch.npnts = ch.U.shape[0]

            elif ch.model == cls.models[4]:
                ch.U, ch.fixedU = cls._get_emo_data(ch.filep)
                ch.npnts = ch.U.shape[0]

            elif ch.model == cls.models[5]:
                ch.U, ch.fixedU = cls._get_mlr_data(ch.filep)
                ch.npnts = ch.U.shape[0]

            elif ch.model == cls.models[6]:
                try:
                    ch.U, ch.fixedU = cls._get_custom_pot_data(ch.filep)
                    ch.npnts = ch.U.shape[0]
                    ch.cfunc = ch.custom_function
                except TypeError as te:
                    print(te)
            else:
                raise SystemExit(f'Invalid potential model {ch.model}')

            if ch.filep not in cls.unique_channels:
                cls.unique_channels[ch.filep] = (ci, ch)
                ch.pointer = (ci, ch)
                ch.isunique = 1
                cls.totp_unique += ch.npnts
            else:
                ch.pointer = cls.unique_channels[ch.filep]
                ch.isunique = 0

            cls.totp += ch.npnts

        cls.ppar = np.array([], dtype=np.float64)
        cls.pfixed = np.array([], dtype=np.int64)
        cls.punits = np.array([], dtype=np.int64)
        for cv in Channel.unique_channels.values():
            cls.ppar = np.append(cls.ppar, cv[1].U)
            # cls.punits = np.append(
            #     cls.punits, np.full(cv[1].U.shape[0], cv[1].yunits)
            # )
            cls.pfixed = np.append(cls.pfixed, cv[1].fixedU)

        # print(cls.unique_channels.keys()) # unique pfiles
        # print(cls.unique_channels.values()) # unique channels

    @classmethod
    def _get_pointwise_data(cls, filep):

        pot = np.loadtxt(filep, skiprows=1)

        R = pot[:, 0] / C_bohr
        U = pot[:, 1]

        fixedU = np.zeros_like(U)
        if pot.shape[1] == 3:
            fixedU = pot[:, 2]

        return R, U, fixedU

    @classmethod
    def _get_morse_data(cls, filep):

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
        U, fixedU = np.zeros(npt), np.zeros(npt, dtype=np.int64)

        U[0] = float(morse_data[mapp[0]][0]) / C_hartree
        U[1] = float(morse_data[mapp[1]][0]) / C_hartree
        U[2] = float(morse_data[mapp[2]][0]) * C_bohr
        U[3] = float(morse_data[mapp[3]][0]) / C_bohr

        if def_fixed:
            for i in range(0, 4):
                fixedU[i] = int(morse_data[mapp[i]][1])

        return U, fixedU

    @classmethod
    def _get_emo_data(cls, filep):

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
            filter(lambda x: x[0].lower().startswith('b'), emo_data.items())
        )

        # TODO: check for key error

        mapp = {0: 'Te', 1: 'De', 2: 'p', 3: 're'}

        U, fixedU = np.zeros(npt), np.zeros(npt, dtype=np.int64)

        U[0] = float(emo_data[mapp[0]][0]) / C_hartree
        U[1] = float(emo_data[mapp[1]][0]) / C_hartree
        U[2] = float(emo_data[mapp[2]][0])
        U[3] = float(emo_data[mapp[3]][0]) / C_bohr

        # all beta coefficients have dimentions 1/distance
        bvalues = list(map(lambda x: float(x[0]) * C_bohr, bparams.values()))
        ni, nb = 4, len(bvalues)

        U[ni:ni+nb] = bvalues

        if def_fixed:
            bfixed = list(map(lambda x: int(x[1]), bparams.values()))

            for i in range(0, ni):
                fixedU[i] = int(emo_data[mapp[i]][1])

            fixedU[ni:ni+nb] = bfixed

        return U, fixedU

    @classmethod
    def _get_mlr_data(cls, filep):
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
            filter(lambda x: x[0].lower().startswith('b'), mlr_data.items())
        )
        cparams = dict(
            filter(lambda x: x[0].lower().startswith('c'), mlr_data.items())
        )
        dparams = dict(
            filter(lambda x: x[0].lower().startswith('d'), mlr_data.items())
        )
        # remove De
        dparams = dict(
            filter(lambda x: x[0].lower() != 'de', dparams.items())
        )
        # remove binf
        bparams = dict(
            filter(lambda x: x[0].lower() != 'binf', bparams.items())
        )
        # TODO: check for key error

        mapp = {
            0: 'Te', 1: 'De', 2: 'p', 3: 'q', 4: 'rref', 5: 're', 6: 'binf'
        }

        U, fixedU = np.zeros(npt), np.zeros(npt, dtype=np.int64)
        U[0] = float(mlr_data[mapp[0]][0]) / C_hartree
        U[1] = float(mlr_data[mapp[1]][0]) / C_hartree
        U[2] = float(mlr_data[mapp[2]][0])
        U[3] = float(mlr_data[mapp[3]][0])
        U[4] = float(mlr_data[mapp[4]][0]) / C_bohr
        U[5] = float(mlr_data[mapp[5]][0]) / C_bohr
        U[6] = float(mlr_data[mapp[6]][0]) * C_bohr

        # TODO: check dimenstions of C and D parameters - they are not correct!
        bvalues = list(
            map(lambda x: float(x[0]) * C_bohr, bparams.values())
        )
        cvalues = list(
            map(lambda x: float(x[0]), cparams.values())
        )
        dvalues = list(
            map(lambda x: float(x[0]), dparams.values())
        )

        ni, nb, nc, nd = 7, len(bvalues), len(cvalues), len(dvalues)
        U[ni:ni+nb] = bvalues
        U[ni+nb:ni+nb+nc] = cvalues
        U[ni+nb+nc:ni+nb+nc+nd] = dvalues

        if def_fixed:
            bfixed = list(map(lambda x: int(x[1]), bparams.values()))
            cfixed = list(map(lambda x: int(x[1]), bparams.values()))
            dfixed = list(map(lambda x: int(x[1]), bparams.values()))

            for i in range(0, ni):
                fixedU[i] = int(mlr_data[mapp[i]][1])

            fixedU[ni:ni+nb] = bfixed
            fixedU[ni+nb:ni+nb+nc] = cfixed
            fixedU[ni+nb+nc:ni+nb+nc+nd] = dfixed

        # channel.npnts = ni + nb + nc + nd
        return U, fixedU

    @classmethod
    def _get_custom_pot_data(cls, filep):
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

        U, fixedU = np.zeros(npt), np.zeros(npt, dtype=np.int64)

        for i, v in enumerate(data.values()):
            U[i] = float(v[0])
            fixedU[i] = int(v[1])

        return U, fixedU

    @classmethod
    def edit_channel_parameters(cls, ypar, channels):

        # the total number of potential points
        cls.tot_npts, onpnts = 0, 0

        # for channel in cls.unique_channels:
        for cv in Channel.unique_channels.values():
            channel = cv[1]
            if channel.model == cls.models[1] or \
               channel.model == cls.models[2]:
                npnts = channel.U.shape[0]
                xpnts = channel.R * C_bohr
                fix = channel.fixedU
                st, en = onpnts, onpnts + npnts
                np.savetxt(
                    channel.filep,
                    np.column_stack([xpnts, ypar[st:en] * C_hartree, fix]),
                    header=str(npnts),
                    comments='',
                    fmt=['%20.12f', '%25.14f', '%6d']
                )
                cls.tot_npts += npnts
                onpnts = en
            elif channel.model == cls.models[3]:
                npnts = channel.U.shape[0]
                st, en = onpnts, onpnts + npnts
                Te, De, a, re = ypar[st:en]

                with open(channel.filep, 'w') as fp:
                    fp.write(
                        f'Te = {Te*C_hartree:>17.8e}'
                        f'{channel.fixedU[0]:>10}\n'
                    )
                    fp.write(
                        f'De = {De*C_hartree:>17.8e}'
                        f'{channel.fixedU[1]:>10}\n'
                    )
                    fp.write(
                        f'a  = {a:>17.8f}{channel.fixedU[2]:>10}\n'
                    )
                    fp.write(
                        f're = {re*C_bohr:>17.8f}'
                        f'{channel.fixedU[3]:>10}\n'
                    )
                cls.tot_npts += npnts
                onpnts = en
            elif channel.model == cls.models[4]:
                pass
            elif channel.model == cls.models[5]:
                pass
            elif channel.model == cls.models[6]:
                npnts = channel.U.shape[0]
                st, en = onpnts, onpnts + npnts
                params = ypar[st:en]

                with open(channel.filep, 'w') as fp:
                    for i in range(0, npnts):
                        fp.write(
                            f'param{i+1} = {params[i]:>17.8e}'
                            f'{channel.fixedU[i]:>10}\n'
                        )
                cls.tot_npts += npnts
                onpnts = en


class Coupling:

    def __init__(self, **kwargs):

        self.interact = kwargs['interact']
        self.coupling = kwargs['coupling']
        self.model = kwargs['model'].lower()
        self.label = kwargs['label']
        self.multiplier = kwargs.get('multiplier') or 1.0
        self.custom_function = kwargs.get('custom_function')

    @classmethod
    def get_coupling_objects(cls):
        return cls.couplings

    @classmethod
    def _define_coupling_models(cls):

        return {
            1: 'pointwise',
            2: 'cspline',
            3: 'constant',
            4: 'custom',
        }

    @classmethod
    def set_coupling_parameters(cls, couplings, cfile='couplings.yml'):

        cls.models = cls._define_coupling_models()
        cls.cpl_file = cfile
        cls.couplings = couplings
        cls.cpar = np.array([], dtype=np.float64)
        cls.cunits = np.array([], dtype=np.float64)
        cls.cfixed = np.array([], dtype=np.int64)
        cls.cregular = np.array([], dtype=np.float64)
        cls.clambda = np.array([], dtype=np.float64)
        cdata = cls._get_couplings_data()

        for cp in couplings:
            params_str = list(map(str.split, cdata[cp.label]))
            params = np.array([list(map(float, i)) for i in params_str])

            # if string convert to tuples
            if isinstance(cp.coupling, str):
                cp.coupling = (cp.coupling,)
                cp.interact = (cp.interact,)
                cp.multiplier = (cp.multiplier,)

            # convert to lowercase
            cp.coupling = tuple(map(str.lower, cp.coupling))

            # set units
            cp.xunits = 1.0 / C_bohr
            cp.yunits = 1.0

            if 'spin-orbit' in cp.coupling:
                cp.yunits = 1.0 / C_hartree

            if all(item.startswith('lambdad') for item in cp.coupling):
                cp.yunits = C_hartree

            # set parameters for pointwise models
            if cp.model == cls.models[1] or cp.model == cls.models[2]:
                cp.xc = params[:, 0] * cp.xunits
                cp.yc = params[:, 1] * cp.yunits
                cls.cpar = np.append(cls.cpar, cp.yc)
                cls.cunits = np.append(cls.cunits, cp.yunits)

                if params.shape[1] >= 3:
                    cp.fixedp = params[:, 2]
                    cls.cfixed = np.append(cls.cfixed, cp.fixedp)
                if params.shape[1] >= 4:
                    cp.regularp = params[:, 3]
                    cp.lambdai = params[:, 4]
                    cls.cregular = np.append(cls.cregular, cp.regularp)
                    cls.clambda = np.append(cls.clambda, cp.lambdai)

            # set parameters for constant or custom models
            if cp.model == cls.models[3] or cp.model == cls.models[4]:
                cp.yc = params[:, 0]

                if params.shape[1] >= 2:
                    cp.fixedp = params[:, 1]
                if params.shape[1] >= 3:
                    cp.regularp = params[:, 2]
                    cp.lambdai = params[:, 3]

            cp.cfunc = cp.custom_function
            cp.npnts = params.shape[0]

    @classmethod
    def _get_couplings_data(cls):

        try:
            with open(cls.cpl_file, 'r') as inps:
                try:
                    return yaml.load(inps)
                except yaml.YAMLError as exc:
                    # TODO: fails with YAML object has no attribute YAMLError
                    raise SystemExit(exc)
        except IOError as e:
            raise SystemExit(e)

    @classmethod
    def edit_coupling_parameters(cls, ypar, couplings):

        cpl_data = cls._get_couplings_data()
        cpar = Channel.tot_npts

        for cp in couplings:
            new_item = []
            if cp.model == cls.models[1] or cp.model == cls.models[2]:
                units = cp.yunits
                for item in cpl_data[cp.label]:
                    sitem = item.split()
                    if len(sitem) >= 3:
                        new_item.append(
                            f'{float(sitem[0]):10.12f}'
                            f'{float(ypar[cpar]/units):24.14f}'
                            f'{int(sitem[2]):7d}'
                            f'{float(sitem[3]):14.2e}'
                            f'{float(sitem[4]):14.1e}'.lstrip()
                        )
                    else:
                        new_item.append(
                            f'{float(sitem[0]):10.12f}'
                            f'{float(ypar[cpar]/units):24.14f}'
                            f'{int(sitem[2]):7d}'
                        )
                    cpar += 1
            else:
                for item in cpl_data[cp.label]:
                    sitem = item.split()
                    if len(sitem) >= 2:
                        new_item.append(
                            f'{float(ypar[cpar]):18.12f}'
                            f'{int(sitem[1]):3d}'
                            f'{float(sitem[2]):10.1}'.lstrip()
                        )
                    else:
                        new_item.append(
                            f'{float(ypar[cpar]):18.12f}'
                            f'{int(sitem[1]):3d}'.lstrip()
                        )
                    cpar += 1
            cpl_data[cp.label] = new_item

        with _open(cls.cpl_file, 'w', encoding='utf8') as stream:
            yaml.dump(cpl_data, stream=stream)
