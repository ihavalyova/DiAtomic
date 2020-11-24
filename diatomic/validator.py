from os import stat as _stat
from os.path import exists as _exists
from logging import error as _error, warning as _warning
import numpy as np

__all__ = ['Validator']


class Validator:

    def __init__(self):
        pass

    @classmethod
    def validate(cls, molecule_data=None, grid=None,
                 channels=None, couplings=None, fitting=None):

        if molecule_data is not None:

            if molecule_data.masses is None and molecule_data.molecule is None:
                _error(
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
                _error('The required parameter "nisotopes" is missing.')
                raise SystemExit()

            if molecule_data.jrange is not None:
                cls._check_jrange(molecule_data.rot_range)
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
            cls._check_channel_data(channels)
        if couplings is not None:
            cls._check_couplings_data(couplings)
        if fitting is not None:
            cls._check_fitting_data(fitting)

    @classmethod
    def _check_type_iterable(cls, iterable):

        try:
            _ = iter(iterable)
        except TypeError as te:
            _error(f'The type of {iterable} raised an error')
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
            _error(f'Masses={masses} raised an error')
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
            _error(f'Parities={parities} raised an error')
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
            _error('The type of "jrange" raised an error')
            raise SystemExit(te)

        try:
            cls._check_negative_values(jrange)
            cls._check_jrange_values(jmin, jmax)
        except ValueError as ve:
            _error(f'jrange data {jmin}, {jmax} raised an error')
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
            _error(
                'Incorrect type for parameters "rrange" and/or "npoints".'
            )
            raise SystemExit(te)

        try:
            cls._check_negative_values([rmax, rmin, ngrid, rbar])  # rbar?
            cls._check_grid_values(rmin, rmax)
        except ValueError as ve:
            _error(
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
    def check_molecule_levels_parameters(cls, ml):

        try:
            cls._check_negative_values([ml.nch])
        except ValueError as ve:
            print(ve)

        if ml.nch == 0:
            raise ValueError('At least one channel should be defined.')

        cls._check_type_iterable(ml.energy_subset_index)
        cls._check_type_iterable(ml.energy_subset_value)

    @classmethod
    def _check_channel_data(cls, channels):

        for ch in channels:
            if ch.model == 'pointwise':
                if not _exists(ch.filep):
                    raise FileNotFoundError(
                        f'The file {ch.filep} does not exist.'
                    )
                if _stat(ch.filep).st_size == 0:
                    raise OSError(f'The file {ch.filep} is empty.')

                pot_data = np.loadtxt(ch.filep, skiprows=1)

                msg_init = 'The pointwise potential file'
                if pot_data.shape[0] < 5:
                    _warning(
                        f'{msg_init} {ch.filep} has less than 5 parameters'
                    )
                if pot_data.shape[0] < 2:
                    _error(
                        f'{msg_init} {ch.filep} has less than 2 parameters'
                    )
                if pot_data.shape[1] == 2:
                    _warning(
                        f'{msg_init} {ch.filep} has 2 columns'
                    )
                if pot_data.shape[1] < 2:
                    _error(
                        f'{msg_init} {ch.filep} has less than 2 columns'
                    )

            if ch.model == 'custom':
                if not hasattr(ch, '__call__'):
                    _error(
                        'Model is set to custum but ' +
                        'custom function is not provided'
                    )

        # TODO: check if model 'custom' is defined and custom_func is missing

    @classmethod
    def _check_couplings_data(cls, couplings):
        pass

    @classmethod
    def _check_fitting_data(cls, fitting):
        pass
