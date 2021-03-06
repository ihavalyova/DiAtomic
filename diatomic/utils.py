from scipy.constants import atomic_mass as _atomic_mass, angstrom as _angstrom
from scipy.constants import physical_constants as _physical_constants
from shutil import copy2 as _copy2
import datetime
import os

__all__ = ['Utils']

# convert from m^-1 to cm^-1
C_hartree = _physical_constants['hartree-inverse meter relationship'][0]/100.0

# convert from m to angstrom
C_bohr = _physical_constants['Bohr radius'][0] / _angstrom

# devide by m_e to convert from amu to au
C_massau = _atomic_mass / _physical_constants['electron mass'][0]


class Utils:

    @classmethod
    def createBackup(cls, ref_file):

        backup_folder = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), '.backup'
        )

        os.makedirs(backup_folder, exist_ok=True)

        fname, fext = os.path.splitext(ref_file)

        name_bkp = '_'.join([fname, 'bkp', cls.getDatetime()]) + fext

        _copy2(ref_file, os.path.join(backup_folder, name_bkp))

    @classmethod
    def getDatetime(cls):

        now = datetime.datetime.now()
        # milli = str(int(round(now.microsecond / 1.0e3)))  # us to ms

        timelist = [
            now.year, now.month, now.day, now.hour,
            now.minute, now.second, now.microsecond
        ]

        return '_'.join(str(t) for t in timelist)

    @classmethod
    def get_current_dir(cls):

        # will not work if the program is used as a module
        # return os.path.abspath(os.path.dirname(__file__))

        return os.getcwd()

    @classmethod
    def get_plot_dir(cls, key):

        plot_dir = 'plotting'
        func_dir = 'interactions'
        resid_dir = 'full_residuals'
        wavefunc_dir = 'wavefunctions'
        res_dir = 'residuals'

        cls.plot_path = os.path.join(cls.get_current_dir(), plot_dir)
        os.makedirs(cls.plot_path, exist_ok=True)
        cls.func_path = os.path.join(cls.plot_path, func_dir)
        cls.resid_path = os.path.join(cls.plot_path, resid_dir)
        cls.wavefunc_path = os.path.join(cls.plot_path, wavefunc_dir)
        cls.res_path = os.path.join(cls.plot_path, res_dir)

        paths = {
            'plot_path': cls.plot_path,
            'func_path': cls.func_path,
            'resid_path': cls.resid_path,
            'wavefunc_path': cls.wavefunc_path,
            'res_path': cls.res_path
        }

        return paths[key]

    @classmethod
    def get_level_output_dir(cls, key):
        pass
