import datetime
import os
import numpy as np
from scipy.constants import atomic_mass as _atomic_mass, angstrom as _angstrom
from scipy.constants import physical_constants as _physical_constants
from shutil import copy2 as _copy2
from math import sqrt as _sqrt

__all__ = ['Utils']


class Utils:

    # convert from m^-1 to cm^-1
    C_hartree = _physical_constants[
        'hartree-inverse meter relationship'][0]/100.0

    # convert from m to angstrom
    C_bohr = _physical_constants['Bohr radius'][0] / _angstrom

    # devide by m_e to convert from amu to au
    C_massau = _atomic_mass / _physical_constants['electron mass'][0]

    c_boltzmank_str = 'Boltzmann constant in inverse meter per kelvin'
    C_boltzmannk = _physical_constants[c_boltzmank_str][0] * 0.01

    @classmethod
    def calculate_stats(cls, yobs, ycal, yunc, is_weighted):

        ndata = yobs.shape[0]
        diff_square = np.square(yobs - ycal)

        # calculate chi square
        if not is_weighted:
            weights = 1.0 / np.square(yunc)
        else:
            weights = 1.0 / (np.square(yunc) + 0.33 * (diff_square))

        chi2 = np.sum(diff_square * weights) / ndata

        # calculate rms
        rms = _sqrt(np.sum(diff_square) / ndata)

        # calculate dimensionless rms
        rmsd = _sqrt(chi2)

        # calculate mean error
        mean_error = np.sum(np.abs(diff_square)) / ndata

        # calculate the average deviation
        avrg = np.mean(yobs)
        avrg_dev = np.sum(np.abs(ycal - avrg)) / ndata

        # calculate standard deviation
        std_dev = _sqrt(np.sum(np.square(ycal - avrg)) / ndata)

        # calculate the median value - the half values are
        # before and the other halfs after it
        median = np.median(ycal)

        stats = {
            'ndata': ndata,
            'chi2': chi2,
            'rms': rms,
            'rmsd': rmsd,
            'mean_error': mean_error,
            'average_dev': avrg_dev,
            'std_dev': std_dev,
            'median': median
        }

        return stats

    @classmethod
    def print_stats(cls, stats):

        out = ''
        init_line = f'  {39*"-"}\n'

        def make_table_row(key, value):
            out_row = ''
            out_row += init_line
            line = ' | {0:^{width1}} | {1:{width2}.{prec}f} |\n'.format(
                key, value, width1=12, width2=22, prec=5)
            out_row += line

            return out_row

        for key in stats.keys():
            out_row = make_table_row(key, stats[key])
            out += out_row

        out += init_line

        return '\n' + out

    @classmethod
    def create_backup_file(cls, ref_file):

        backup_folder = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), '.backup')

        os.makedirs(backup_folder, exist_ok=True)

        fname, fext = os.path.splitext(ref_file)

        name_bkp = '_'.join([fname, 'bkp', cls.getDatetime()]) + fext

        _copy2(ref_file, os.path.join(backup_folder, name_bkp))

    @classmethod
    def get_date_time(cls):

        now = datetime.datetime.now()
        # milli = str(int(round(now.microsecond / 1.0e3)))  # us to ms

        timelist = [now.year, now.month, now.day, now.hour,
                    now.minute, now.second, now.microseconds]

        return '_'.join(str(t) for t in timelist)

    @classmethod
    def get_current_directory(cls):

        # will not work if the program is used as a module
        return os.path.abspath(os.path.dirname(__file__))

        # return os.getcwd()

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
    def print_detailed_output(cls, H, file_name=None):

        file_name = file_name or 'data_info.dat'
        nsymb = 50
        s = ' '

        mass_str = ''
        for n, nisotope in enumerate(H.niso):
            mass = H.masses[n]
            mass_str += (
                f'{s:>4}Isotopologue {nisotope}: '
                f'Reduced Mass = {mass:>15.9} au = '
                f'{mass / Utils.C_massau:>15.9} amu\n')

        jqnum_str = ' '.join(map(str, H.jqnumbers))
        pars_str = ', '.join(map(lambda x: 'e' if x else 'f', H.pars))
        shift_enr = True if H.refj is not None else False
        iso_str = ' '.join(map(str, H.molecule))

        output_mol_data = (
            f'{nsymb*"#"} Molecule Data {nsymb*"#"}\n\n'
            f'{s:>4}Chemical Symbol of the Molecule: {H.mol_name}\n\n'
            f'{s:>4}Number of Defined Isotopologues = '
            f'{len(H.masses):>4d}\n\n'
            f'{s:>4}Isotopologues = {iso_str}\n\n'
            f'{s:>4}Number of Used Isotopologues = {len(H.niso):>7d}\n\n'
            f'{mass_str}\n\n'
            f'{s:>4}Number of Rotational Quantum Numbers = '
            f'{H.jqnumbers.shape[0]}\n\n'
            f'{s:>4}Rotational Quantum Numbers:\n'
            f'{s:>4}{jqnum_str}\n\n'
            f'{s:>4}Symmetry Levels = '
            f'{pars_str}\n\n'
            f'{s:>4}Shift Energies = {shift_enr}\n'
            f'{s:>4}Shift Energies by Level J = {H.refj}\n\n')

        grid_points_str = ''
        for point in H.rgrid:
            grid_points_str += (
                f'{s:<4}{point*Utils.C_bohr:>20.8f};'
                f'{point:>20.8f}\n')

        output_grid = (
            f'{nsymb*"#"} Grid {nsymb*"#"}\n\n'
            f'{s:>4}Method of Solution: {H.solver}\n\n'
            f'{s:>4}Number of Grid Points = {H.ngrid:<5d}\n\n'
            f'{s:>4}Rmin = {H.rmin*Utils.C_bohr:>12.10f} '
            f'Angstrom = {H.rmin:>12.10f} Bohr\n'
            f'{s:>4}Rmax = {H.rmax*Utils.C_bohr:>12.10f} '
            f'Angstrom = {H.rmax:>12.10f} Bohr\n\n'
            f'{s:>4}Hamiltonian Matrix Size = '
            f'{s:>4}{H.nch*H.ngrid} x {H.nch*H.ngrid}\n'
            f'{s:>4}Number of Computed Eigenvalues = '
            f'{H.ecount:>8d}\n'
            f'{s:>4}Number of Selected Eigenvalues = '
            f'{H.out_data.shape[0]:>8d}\n\n'
            f'{s:>4}Grid Points (Angstrom; Bohr):\n\n'
            f'{grid_points_str}\n\n')

        channels_str = ''
        for ic, ch in enumerate(H.channels):
            eq_pos = np.argmin(ch.upoints)
            channels_str += (
                f'\n{s:>9} {ic+1}.\n'
                f'{s:<13}Model: {ch.model}\n'
                f'{s:<13}File: {ch.filep}\n'
                f'{s:<13}Lambda: {ch.nlambda}\n'
                f'{s:<13}Sigma: {ch.nsigma}\n'
                f'{s:<13}Multiplicity: {ch.mult}\n'
                f'{s:<13}Rot correction: {ch.rot_correction}\n'
                f'{s:<13}Equilibrium distance point = {eq_pos+1}\n'
                f'{s:<13}Equilibrium distance: '
                f'Rmin = {ch.rpoints[eq_pos]*Utils.C_bohr}, '
                f'Umin = {ch.upoints[eq_pos]*Utils.C_hartree}\n'
                f'{s:<13}Number of parameters: {ch.npnts}\n'
                f'{s:<13}Parameters (Angstrom/cm-1; Bohr/Hartree):\n\n')

            if ch.model == 'pointwise':
                for i in range(0, len(ch.rpoints)):
                    channels_str += (
                        f'{s:<4}{ch.rpoints[i]*Utils.C_bohr:>20.8f}'
                        f'{ch.upoints[i]*Utils.C_hartree:>20.8f};'
                        f'{ch.rpoints[i]:>20.8f}'
                        f'{ch.upoints[i]:>20.8f}\n')

        output_channels = (
            f'{nsymb*"#"} Channels Data {nsymb*"#"}\n\n'
            f'{s:>4}Number of Channels = {H.nch}\n'
            f'{channels_str}')

        ugrid_cols = np.hstack((
            H.rgrid[:, np.newaxis] * Utils.C_bohr,
            H.ugrid.reshape(H.nch, H.ngrid).T * Utils.C_hartree))

        output_channels_funcs = (
            f'\n{nsymb*"#"} Channel Functions on Grid {nsymb*"#"}\n\n'
            f'{ugrid_cols}\n')

        couplings_str = ''
        for ic, cp in enumerate(H.couplings):
            couplings_str += (
                f'{s:>9} {ic+1}.\n'
                f'{s:<13}Type: {cp.model}\n'
                f'{s:<13}Channels: {cp.interact}\n'
                f'{s:<13}Coupling: {cp.coupling}\n'
                f'{s:<13}Label: {cp.label}\n'
                f'{s:<13}Multiplier: {cp.multiplier}\n'
                f'{s:<13}Parameters:\n')

            if cp.model == 'pointwise':
                for i in range(0, len(cp.xc)):
                    couplings_str += (
                        f'{cp.xc[i]:>29.14f}'
                        f'{cp.yc[i]:>25.14f}\n')

        output_couplings = (
            f'{nsymb*"#"} Couplings Data {nsymb*"#"}\n\n'
            f'{s:>4}Number of Couplings = {H.ncp}\n\n'
            f'{couplings_str}')

        fgrid_cols = np.hstack((
            H.rgrid[:, np.newaxis] * Utils.C_bohr,
            H.fgrid.reshape(H.ncp, H.ngrid).T))

        output_couplings_funcs = (
            f'{nsymb*"#"} Coupling Functions on Grid {nsymb*"#"}\n\n'
            f'{fgrid_cols}\n')

        output_exp_energies = ''
        if H.exp_data is not None:
            output_exp_energies = (
                f'{nsymb*"#"} Experimental data {nsymb*"#"}\n\n'
                f'{s:>4}File with Experimental Data = {H.exp_file}\n\n'
                f'{s:>4}Markers: \n'
                f'{s:>4}Number of used experimental data = '
                f'{H.exp_data.shape[0]}\n\n')

        with open(file_name, 'w') as outf:
            outf.write(output_mol_data)
            outf.write(output_grid)
            outf.write(output_channels)
            outf.write(output_channels_funcs)
            outf.write(output_couplings)
            outf.write(output_couplings_funcs)
            outf.write(output_exp_energies)
