import io
import os
import re
import logging
import numpy as np
from more_itertools import unique_everseen
from collections import OrderedDict
from atomic_database import AtomicDatabase
from constants import Const

# to check if C loader is present on current machine
# /path/to/python -c "import ruamel.yaml;
# print(ruamel.yaml.__with_libyaml__)"

from ruamel.yaml import YAML
yaml = YAML(typ='rt', pure=False)


class MoleculeData:

    def __init__(self):

        self.imasses = None
        self.molecule_names = None
        self.pars = (0, 1)
        self.niso = None
        self.rots = None
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

        self.molecule_names = molecule_names

        for molecule_name in self.molecule_names:
            rmass = self._calculate_molecule_reduced_mass(molecule_name)
            rmass *= Const.massau
            self.reduced_masses.append(rmass)

    @property
    def masses(self):

        return self.imasses

    @masses.setter
    def masses(self, imasses):

        self.imasses = [m * Const.massau for m in imasses]

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
    def jrange(self, rots):

        self.rots = rots

        self.jqnumbers = np.unique(np.hstack(
            (
                self.jqnumbers,
                np.arange(rots[0], rots[1]+1.0, dtype=np.float64)
            )
        ))

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

        self.refE = refE / Const.hartree

    def _arange_jqnumbers(self):

        self.jqnumbers = np.delete(
            self.jqnumbers, np.argwhere(self.jqnumbers == self.refj)
        )

        self.jqnumbers = np.insert(self.jqnumbers, 0, self.refj)

    def _calculate_molecule_reduced_mass(self, molecule, database=None):

        found = re.search(r'^(\d+)(\D+)(\d+)(\D+)\b', molecule.strip())

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

        pattern = re.compile(
            fr'Atomic\s+Symbol\s+=\s+{symbol}[\r\n]'
            fr'Mass\s+Number\s+=\s+{mnumber}'
            fr'[\r\n]Relative\s+Atomic\s+Mass\s+=\s+\d+\.\d+'
        )

        atom_data = re.findall(pattern, atomic_db)

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
            return self._read_experimental_data(markers, average)
        except Exception as e:
            raise SystemExit(e)

    def _read_experimental_data(self, markers, average):

        ndata = np.genfromtxt(
            self.exp_file, max_rows=1, comments='#'
        )

        self.exp_data = np.genfromtxt(
            self.exp_file, skip_header=1, max_rows=int(ndata), comments='#'
        )

        # filter by marker
        if markers is not None:
            marker_mask = np.in1d(
                self.exp_data[:, 6], np.fromiter(markers, dtype=np.int64)
            )
            self.exp_data = self.exp_data[marker_mask]

        # filter by parity
        if len(self.pars) == 1:
            parity_mask = np.in1d(
                self.exp_data[:, 4], np.fromiter(self.pars, dtype=np.int64)
            )
            self.exp_data = self.exp_data[parity_mask]

        self.jqnumbers = np.intersect1d(
            self.jqnumbers, np.unique(self.exp_data[:, 2])
        )

        if self.refj is not None:
            self._arange_jqnumbers()

        # TODO: change parity above in the same way as j

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

        fname, fext = os.path.splitext(self.exp_file)
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
        self.rot_correction = kwargs.get('rot_correction') or 0.0
        self.custom_function = kwargs.get('custom_function')

    @classmethod
    def get_unique_channels_parameters(cls, channels):

        pars = np.array([], dtype=np.float64)
        fixed = np.array([], dtype=np.int64)

        cls.unique_pfiles = OrderedDict()
        cls.unique_channels = OrderedDict()

        for channel in channels:
            if channel.filep not in cls.unique_pfiles:

                if channel.model == cls.models[1] or \
                  channel.model == cls.models[2]:
                    _, U, fixedU = cls._get_pointwise_data(channel.filep)
                elif channel.model == cls.models[3]:
                    U, fixedU = cls._get_morse_data(channel.filep)
                elif channel.model == cls.models[4]:
                    U, fixedU = cls._get_emo_data(channel.filep)
                elif channel.model == cls.models[5]:
                    U, fixedU = cls._get_mlr_data(channel.filep)
                elif channel.model == cls.models[6]:
                    U, fixedU = cls._get_custom_pot_data(channel.filep)
                elif channel.model == cls.models[7]:
                    pass

                pars = np.append(pars, U)
                fixed = np.append(fixed, fixedU)

                cls.unique_pfiles[channel.filep] = channel.filep
                cls.unique_channels[channel] = channel

        cls.unique_pfiles = list(cls.unique_pfiles)
        cls.unique_channels = list(cls.unique_channels)

        return pars, fixed

    @classmethod
    def _get_channel_models(cls):

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

        cls.models = cls._get_channel_models()

        for channel in channels:

            if channel.model == cls.models[1] or \
               channel.model == cls.models[2]:

                channel.R, channel.U, channel.fixedU = \
                    cls._get_pointwise_data(channel.filep)

                channel.npnts = channel.U.shape[0]

            elif channel.model == cls.models[3]:

                channel.U, channel.fixedU = \
                    cls._get_morse_data(channel.filep)
                channel.npnts = channel.U.shape[0]

            elif channel.model == cls.models[4]:

                channel.U, channel.fixedU = \
                    cls._get_emo_data(channel.filep)
                channel.npnts = channel.U.shape[0]

            elif channel.model == cls.models[5]:

                channel.U, channel.fixedU = \
                    cls._get_mlr_data(channel.filep)
                channel.npnts = channel.U.shape[0]

            elif channel.model == cls.models[6]:

                try:
                    channel.U, channel.fixedU = \
                        cls._get_custom_pot_data(channel.filep)

                    channel.npnts = channel.U.shape[0]
                    channel.cfunc = channel.custom_function
                except TypeError as te:
                    print(te)
            else:
                raise SystemExit(f'Invalid potential model {channel.model}')

        # TODO: remove this
        cls.pnames = list(
            unique_everseen([ch.filep for i, ch in enumerate(channels)])
        )

    @classmethod
    def _get_pointwise_data(cls, filep):

        pot = np.loadtxt(filep, skiprows=1)
        R = pot[:, 0] / Const.bohr
        U = pot[:, 1] / Const.hartree

        # set default third column if it's not provided
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

        U[0] = float(morse_data[mapp[0]][0]) / Const.hartree
        U[1] = float(morse_data[mapp[1]][0]) / Const.hartree
        U[2] = float(morse_data[mapp[2]][0]) * Const.bohr
        U[3] = float(morse_data[mapp[3]][0]) / Const.bohr

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

        U[0] = float(emo_data[mapp[0]][0]) / Const.hartree
        U[1] = float(emo_data[mapp[1]][0]) / Const.hartree
        U[2] = float(emo_data[mapp[2]][0])
        U[3] = float(emo_data[mapp[3]][0]) / Const.bohr

        # all beta coefficients have dimentions 1/distance
        bvalues = list(
            map(lambda x: float(x[0]) * Const.bohr, bparams.values())
        )
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

        U[0] = float(mlr_data[mapp[0]][0]) / Const.hartree
        U[1] = float(mlr_data[mapp[1]][0]) / Const.hartree
        U[2] = float(mlr_data[mapp[2]][0])
        U[3] = float(mlr_data[mapp[3]][0])
        U[4] = float(mlr_data[mapp[4]][0]) / Const.bohr
        U[5] = float(mlr_data[mapp[5]][0]) / Const.bohr
        U[6] = float(mlr_data[mapp[6]][0]) * Const.bohr

        # TODO: check dimenstions of C and D parameters - they are not correct!

        bvalues = list(
            map(lambda x: float(x[0]) * Const.bohr, bparams.values())
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

        for channel in cls.unique_channels:

            if channel.model == cls.models[1] or \
              channel.model == cls.models[2]:

                npnts = channel.U.shape[0]
                xpnts = channel.R * Const.bohr
                fix = channel.fixedU
                st, en = onpnts, onpnts + npnts

                np.savetxt(
                    channel.filep,
                    np.column_stack([xpnts, ypar[st:en] * Const.hartree, fix]),
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
                        f'Te = {Te*Const.hartree:>17.8e}'
                        f'{channel.fixedU[0]:>10}\n'
                    )
                    fp.write(
                        f'De = {De*Const.hartree:>17.8e}'
                        f'{channel.fixedU[1]:>10}\n'
                    )
                    fp.write(
                        f'a  = {a:>17.8f}{channel.fixedU[2]:>10}\n'
                    )
                    fp.write(
                        f're = {re*Const.bohr:>17.8f}{channel.fixedU[3]:>10}\n'
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
    def _get_coupling_models(cls):

        return {
            1: 'pointwise',
            2: 'cspline',
            3: 'constant',
            4: 'custom',
        }

    @classmethod
    def set_coupling_parameters(cls, cpl_file, couplings):

        cls.models = cls._get_coupling_models()
        cls.cpl_file = cpl_file
        coupling_data = cls._get_couplings_data()

        # convert integer keys to string
        coupling_data = {str(k): v for k, v in coupling_data.items()}

        for coupling in couplings:

            params = np.array(
                list(map(str.split, coupling_data[coupling.label]))
            )

            if coupling.model == cls.models[1] or \
               coupling.model == cls.models[2]:
                coupling.xc = np.fromiter(
                    map(float, params[:, 0]), dtype=np.float64
                ) / Const.bohr
                coupling.yc = np.fromiter(
                    map(float, params[:, 1]), dtype=np.float64
                )
                if params.shape[1] >= 3:
                    coupling.fixedp = np.fromiter(
                        map(int, params[:, 2]), dtype=np.int64
                    )
                if params.shape[1] >= 4:
                    coupling.regularp = np.fromiter(
                        map(float, params[:, 3]), dtype=np.float64
                    )
                    coupling.lambdai = np.fromiter(
                        map(float, params[:, 4]), dtype=np.float64
                    )

            if coupling.model == cls.models[3] or \
               coupling.model == cls.models[4]:

                coupling.yc = np.fromiter(
                    map(float, params[:, 0]), dtype=np.float64
                )
                if params.shape[1] >= 2:
                    coupling.fixedp = np.fromiter(
                        map(int, params[:, 1]), dtype=np.int64
                    )
                if params.shape[1] >= 3:
                    coupling.regularp = np.fromiter(
                        map(float, params[:, 2]), dtype=np.float64
                    )
                    coupling.lambdai = np.fromiter(
                        map(float, params[:, 3]), dtype=np.float64
                    )

                coupling.cfunc = coupling.custom_function

            coupling.npnts = params.shape[0]

    @classmethod
    def _get_couplings_data(cls):

        try:
            with open(cls.cpl_file, 'r') as inps:
                try:
                    return yaml.load(inps)
                except yaml.YAMLError as exc:
                    raise SystemExit(exc)
        except IOError as e:
            raise SystemExit(e)

    @classmethod
    def get_coupling_parameters(cls, couplings):

        coupling_params = np.array([], dtype=np.float64)
        coupling_fixed = np.array([], dtype=np.int64)
        coupling_regular = np.array([], dtype=np.float64)
        coupling_lambda = np.array([], dtype=np.float64)

        coupling_data = cls._get_couplings_data()

        # convert integer keys to string
        coupling_data = {str(k): v for k, v in coupling_data.items()}

        for coupling in couplings:
            params = np.array(
                list(map(str.split, coupling_data[coupling.label]))
            )

            if coupling.model == cls.models[1] or \
               coupling.model == cls.models[2]:
                coupling_params = np.append(
                    coupling_params,
                    np.fromiter(
                        map(float, params[:, 1]), dtype=np.float64
                    )
                )
                coupling_fixed = np.append(
                    coupling_fixed,
                    np.fromiter(
                        map(int, params[:, 2]), dtype=np.int64
                    )
                )
                if params.shape[1] >= 3:
                    coupling_regular = np.append(
                        coupling_regular,
                        np.fromiter(
                            map(float, params[:, 3]), dtype=np.float64
                        )
                    )
                    coupling_lambda = np.append(
                        coupling_lambda,
                        np.fromiter(
                            map(float, params[:, 4]), dtype=np.float64
                        )
                    )
            else:
                coupling_params = np.append(
                    coupling_params,
                    np.fromiter(
                        map(float, params[:, 0]), dtype=np.float64
                    )
                )
                coupling_fixed = np.append(
                    coupling_fixed,
                    np.fromiter(
                        map(int, params[:, 1]), dtype=np.int64
                    )
                )
                if params.shape[1] >= 2:
                    coupling_regular = np.append(
                        coupling_regular,
                        np.fromiter(
                            map(float, params[:, 2]), dtype=np.float64
                        )
                    )
                    coupling_lambda = np.append(
                        coupling_lambda,
                        np.fromiter(
                            map(float, params[:, 3]), dtype=np.float64
                        )
                    )

        return coupling_params, coupling_fixed, \
            coupling_regular, coupling_lambda

    @classmethod
    def edit_coupling_parameters(cls, ypar, couplings):

        cpl_data = cls._get_couplings_data()

        # convert integer keys to string
        cpl_data = {str(k): v for k, v in cpl_data.items()}

        cpar = Channel.tot_npts

        for coupling in couplings:
            new_item = []

            if coupling.model == cls.models[1] or \
               coupling.model == cls.models[2]:

                for item in cpl_data[coupling.label]:
                    sitem = item.split()

                    if len(sitem) >= 3:
                        new_item.append(
                            f'{float(sitem[0]):10.12f}'
                            f'{float(ypar[cpar]):24.14f}'
                            f'{int(sitem[2]):7d}'
                            f'{float(sitem[3]):14.2e}'
                            f'{float(sitem[4]):14.1e}'.lstrip()
                        )
                    else:
                        new_item.append(
                            f'{float(sitem[0]):10.12f}'
                            f'{float(ypar[cpar]):24.14f}'
                            f'{int(sitem[2]):7d}'
                        )

                    cpar += 1
            else:
                for item in cpl_data[coupling.label]:
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

            cpl_data[coupling.label] = new_item

        with io.open(cls.cpl_file, 'w', encoding='utf8') as stream:
            yaml.dump(cpl_data, stream=stream)


class Validator:

    def __init__(self):
        pass

    @classmethod
    def validate(cls, molecule_data=None, grid=None,
                 channels=None, couplings=None, fitting=None):

        if molecule_data is not None:

            if molecule_data.masses is None and molecule_data.molecule is None:
                logging.error(
                    f'At least one of the required parameters '
                    f'{molecule_data.massses} or {molecule_data.molecule}'
                    f'should be set.'
                )
                raise SystemExit()

            if molecule_data.masses is not None:
                cls._check_masses(molecule_data.masses)

            if molecule_data.molecule is not None:
                cls._check_molecule(molecule_data.molecule)

            if molecule_data.nisotopes is None:
                logging.error('The required parameter "nisotopes" is missing.')
                raise SystemExit()

            if molecule_data.jrange is not None:
                cls._check_jrange(molecule_data.rots)
                cls._check_jqnumbers(molecule_data.jqnumbers)

            if molecule_data.jvalues is not None:
                cls._check_jvalues(molecule_data.rot_values)
                cls._check_jqnumbers(molecule_data.jqnumbers)

            if molecule_data.referencej is not None:
                cls._check_referencej(molecule_data.referencej)

            if molecule_data.parities is not None:
                cls._check_parity(molecule_data.parities)
        if grid is not None:
            # TODO: split this to several methods
            cls._check_grid_data(
                grid.rgrid, grid.npoints, grid.alpha, grid.rbar
            )
        if channels is not None:
            pass
        if couplings is not None:
            pass
        if fitting is not None:
            pass

    @classmethod
    def _check_type_iterable(cls, iterable):

        try:
            _ = iter(iterable)
        except TypeError as te:
            logging.error(f'The type of {iterable} raised an error')
            raise SystemExit(te)
        else:
            if isinstance(iterable, str):
                raise SystemExit(
                    f'The type of {iterable} must not be a string.'
                )

    @classmethod
    def _check_masses(cls, masses):

        cls._check_type_iterable(masses)

        try:
            cls._check_negative_values(masses)
        except ValueError as ve:
            logging.error(f'Masses={masses} raised an error')
            raise SystemExit(ve)

    @classmethod
    def _check_molecule(cls, molecule):

        cls._check_type_iterable(molecule)

    @classmethod
    def _check_parity(cls, parities):

        try:
            cls._check_negative_values(parities)
            cls._check_parity_values(parities)
            cls._check_parity_length(parities)
        except (ValueError, IndexError) as e:
            logging.error(f'Parities={parities} raised an error')
            raise SystemExit(e)

    @classmethod
    def _check_parity_length(cls, parities):

        if len(parities) > 2:
            raise IndexError(
                f'Only One or two parity={parities} values are allowed.'
            )

    @classmethod
    def _check_parity_values(cls, parities):

        if not set(parities).issubset(set([0, 1])):
            raise ValueError(
                f'0 and/or 1 are the only possible parity={parities} values.'
            )

    @classmethod
    def _check_jrange(cls, jrange):

        cls._check_type_iterable(jrange)

        try:
            cls._check_jrange_length(jrange)
        except IndexError as ie:
            print(ie)

        jmin, jmax = jrange[0], jrange[1]

        try:
            for x in [jmin, jmax]:
                cls._check_type_int_or_float(x)
        except TypeError as te:
            logging.error('The type of "jrange" raised an error')
            raise SystemExit(te)

        try:
            cls._check_negative_values(jrange)
            cls._check_jrange_values(jmin, jmax)
        except ValueError as ve:
            logging.error(f'jrange data {jmin}, {jmax} raised an error')
            raise SystemExit(ve)

    @classmethod
    def _check_jrange_values(self, jmin, jmax):

        if jmin >= jmax:
            raise ValueError(
                f'jmax={jmax} should be greater than jmin={jmin}.'
            )

    @classmethod
    def _check_jrange_length(self, jrange):

        if len(jrange) != 2:
            raise IndexError(
                f'Two values should be specfied for jrange={jrange} parameter'
            )

    @classmethod
    def _check_jvalues(cls, jvalues):
        pass

    @classmethod
    def _check_jqnumbers(cls, jqnumbers):

        if len(jqnumbers) == 0:
            raise SystemExit(
                'Error: The properties jrange/jvalues are not set'
            )

    @classmethod
    def _check_referencej(cls, refj):
        cls._check_type_int_or_float(refj)

    @classmethod
    def _check_grid_data(cls, rrange, ngrid, alpha, rbar):

        cls._check_type_iterable(rrange)

        try:
            cls._check_rrange_len(rrange)
        except IndexError as ie:
            print(ie)

        rmin, rmax = rrange[0], rrange[1]

        try:
            for x in [ngrid, rmin, rmax, rbar]:
                cls._check_type_int_or_float(x)
        except TypeError as te:
            logging.error(
                'Incorrect type for parameters "rrange" and/or "npoints".'
            )
            raise SystemExit(te)

        try:
            cls._check_negative_values([rmax, rmin, ngrid, rbar])  # rbar?
            cls._check_grid_values(rmin, rmax)
        except ValueError as ve:
            logging.error(
                f'Incorrect values for {rmin}//{rmax}//{ngrid}.'
            )
            raise SystemExit(ve)

    @classmethod
    def _check_rrange_len(cls, rrange):

        if len(rrange) != 2:
            raise IndexError(
                f'Two values should be specified for rrange={rrange}.'
            )

    @classmethod
    def _check_grid_values(cls, rmin, rmax):

        if rmin >= rmax:
            raise ValueError(
                f'rmax={rmax} should be greater than rmin={rmin}.'
            )

    @classmethod
    def _check_type_int_or_float(cls, param):

        if not isinstance(param, (int, float)):
            raise TypeError(
                f'Parameter value {param} should be an integer or float number'
            )

    @classmethod
    def _check_negative_values(cls, params):

        if any(x < 0 for x in params):
            raise ValueError(
                'Negative parameter values are not allowed.'
            )

    @classmethod
    def _check_channel_data(cls, ml):

        try:
            cls._check_negative_values([ml.nch])
        except ValueError as ve:
            print(ve)

        if ml.nch == 0:
            raise ValueError('At least one channel should be defined.')

        cls._check_type_iterable(ml.energy_subset_index)

        cls._check_type_iterable(ml.energy_subset_value)

    # TODO: check if model 'custom' is defined and custom_func is missing
