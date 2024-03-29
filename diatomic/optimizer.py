import io
from os.path import join as _join
import scipy as sp
import numpy as np

from .utils import Utils as _utils

try:
    from iminuit import Minuit
except ModuleNotFoundError:
    pass
    # print("'iminuit' module is not installed!\n")

__all__ = ['Optimizer']


class Optimizer:
    """Optimize of the radial parameters based on the  provided observed data for
       the energy levels, transition frequencies or transition intensities.
    """

    def __init__(self, H, S=None, output_progress=False, is_weighted=False):

        self.H = H
        self.S = S
        self.progress = output_progress
        self.is_weighted = is_weighted
        self.progress_str = io.StringIO()
        self.sing_values = None

    def _calculate_levels(self, ypar):

        out_levels = self.H.wrapper_calculate_diatomic_levels(ypar)
        ycal = out_levels[:, 7]  # / C_hartree
        yobs = out_levels[:, 8]  # / C_hartree
        ydel = -out_levels[:, 9]  # / C_hartree
        yvar = out_levels[:, 10]  # / C_hartree

        stats = _utils.calculate_stats(yobs, ycal, yvar, self.is_weighted)

        return ycal, yobs, ydel, yvar, stats

    def _calculate_wavenumbers(self, ypar):

        if self.S is None:
            raise ValueError('The "Spectrum" object should be passed as '
                             'an argument to the __init__() method')

        out_wavens = self.S.wrapper_calculate_wavenumbers(ypar)
        ycal = out_wavens[:, 16]  # / C_hartree
        yobs = out_wavens[:, 17]  # / C_hartree
        ydel = out_wavens[:, 18]  # / C_hartree
        yvar = out_wavens[:, 19]  # / C_hartree

        stats = _utils.calculate_stats(yobs, ycal, yvar, self.is_weighted)

        return ycal, yobs, ydel, yvar, stats

    def _calculate_intensity(self, ypar):

        if self.S is None:
            raise ValueError('The "Spectrum" object should be passed as '
                             'an argument to the __init__() method')

        acoefs = self.S.wrapper_calculate_Einstein_coefficients(ypar)
        ycal = acoefs[:, 23]
        yobs = acoefs[:, 20]
        ydel = yobs-ycal
        yvar = acoefs[:, 21]

        stats = _utils.calculate_stats(yobs, ycal, yvar, self.is_weighted)

        print('stats intensity', stats)

        return ycal, yobs, ydel, yvar, stats

    def optimize(self, ypar, yfixed, calculate_data, niter, derivative,
                 tol, lapack_driver, step_size):
        """Implementation of the SVD algorithm used for nonliner least squares fitting.
        Solves the linear system Ax = b using SVD. Starting with a trial set of
        parameters, the system Ax=b is solved on each iteration for the
        corrections x to the current set of parameters x_0, yielding and a new
        set of parameters x_1 so that: (x_1 = x_0 + x). The matrix A contains
        the derivatives of the data wrt the parameters and b is the rhs vector.

        Args:
            ypar (array): The parameters which to be optimized
            yfixed (array): Whether the parameter is fixed or free
            calculate_data ([type]): [description]
            niter (int): The number of iterations
            derivative (str): Whether to use numerical or analytical derivative
            tol (float): The tolerance parameter of SVD: determines which
                linear combinations of parameters to be discarded by SVD
            lapack_driver (str): [description]
            step_size (float): Used to change the parameters during
                the computation of derivatives
        Returns:
            array: The optimized parameters
        """

        chisq_best, rms_best = 0.0, 0.0
        eps = 1.0e-3

        print(f'\nTotal number of parameters = {ypar.shape[0]}'
              f'\nTotal number of parameters to be optimized = '
              f'{yfixed[yfixed == 1].shape[0]}\n')

        for itr in range(1, niter+1):

            # get the energies with the initial/updated parameters from
            # the first/previous iteration
            ycal, yobs, ydel, yvar, stats = calculate_data(ypar)
            chisq_best, rms_best = stats['chi2'], stats['rms']

            print(self.print_info(itr, tol, stats))

            # build the matrix of the derivatives w.r.t. parameters
            if derivative == 'n':
                dydp = self._calculate_numerical_derivatives(
                    ypar, yfixed, calculate_data, ycal, step_size=step_size)
                # np.savetxt('dydp_numerical.dat', dydp, fmt='%.6e')
            else:
                pass

            # calculate the data matrix A from derivative matrix
            data_matrix = (dydp.T / yvar).T

            # calculate the rhs vector b
            rhs_vector = ydel / yvar

            # the propossed corrections to the parameters is the solution
            x, rank, sing_values = self._find_corrections(
                data_matrix, rhs_vector, tol, lapack_driver)

            # calculate new energies using the proposed corrections

            *_, stats_try = calculate_data(ypar+x)
            chi2_try, rms_try = stats_try['chi2'], stats_try['rms']

            # try to change the corrections
            if chi2_try + eps < chisq_best:
                chisq_best, rms_best = chi2_try, rms_try
                ypar = ypar + x
                is_failed = False
            else:  # chi2try >= chi2best
                step = 10.0
                ypar_try = np.copy(ypar)
                x_try = np.copy(x)
                x_try = x_try / step
                is_failed = True

                while True:
                    ypar_try = ypar_try + x_try

                    *_, stats_new = calculate_data(ypar_try)
                    chi2_new, rms_new = stats_new['chi2'], stats_new['rms']

                    if chi2_new + eps < chisq_best:
                        chisq_best, rms_best = chi2_new, rms_new
                        ypar = ypar_try
                        is_failed = False
                        print(f'Smaller step => Chi square = '
                              f'{chisq_best:.6f} | RMS = {rms_best:.6f}')
                    else:  # chi2new >= chi2best:
                        break

            # on odd iteration try to change tol if nothing has been improved
            if itr > 1 and itr % 2 == 1 and is_failed:
                tol /= 5.0
                print(f'\nTol changed = {tol:.1e}')

            # write detailed fit progress
            if self.progress:
                self._write_fit_progress(rank, sing_values, tol, ypar, x, yfixed)
                print(self.progress_str.getvalue())

        self.sing_values = sing_values

        return ypar

    def optimize_levels_by_svd(self, niter=1, derivative='n', tol=0.1,
                               lapack_driver='gelsd', step_size=1.0e-4):
        """This procedure applies the SVD algorithm to optimize the parameters
        using an input information about the energy levels

        Args:
            niter (int, optional): The number of iterations. Defaults to 1.
            derivative (str, optional): Type of derivative. Defaults to 'n'.
            tol (float, optional): The tolerance value. Defaults to 0.1.
            lapack_driver (str, optional): Which lapack procedure to use.
                Defaults to 'gelsd'.
            step_size (float, optional): The value used to change the
                parameters. Defaults to 1.0e-4.

        Returns:
            array: The optimized parameters
        """

        ypar_final = self.optimize(self.H.params, self.H.fixed, self._calculate_levels,
                                   niter, derivative, tol, lapack_driver, step_size)
        self._store_output_params()

        return ypar_final

    def _store_output_params(self):

        prefix = 'fit_'
        for opr in self.H.ops:
            fname = prefix + self.H.fname_data_params
            if opr.model in ['pointwise', 'cspline']:
                # TODO: update ypoints with the fitted params
                self._write_pointwise_data(fname, opr.xpoints, opr.ypoints, opr.fixed)
            elif opr.model == 'morse':
                pass
            elif opr.model == 'emo':
                pass
            elif opr.model == 'mlr':
                pass
            elif opr.model == 'custom':
                pass

    def _write_pointwise_data(self, fname, xpoints, ypoints, fixed):

        data = np.column_stack((xpoints, ypoints, fixed))
        fmt = ['%20.10f', '%25.10f', '%6d']
        np.savetxt(fname, data, fmt=fmt)

    def _write_morse_data(self, fname, ypoints, fixed):

        with open(fname, 'w', encoding='utf8') as fs:
            fs.write(f'Te = {ypoints[0]:>17.8e}{fixed[0]:>9}\n')
            fs.write(f'De = {ypoints[1]:>17.8e}{fixed[1]:>9}\n')
            fs.write(f'a  = {ypoints[2]:>17.8f}{fixed[2]:>9}\n')
            fs.write(f're = {ypoints[3]:>17.8f}{fixed[3]:>9}\n')

    def optimize_wavenumbers_by_svd(self, niter=1, derivative='n', tol=0.1,
                                    lapack_driver='gelsd', step_size=1.0e-4):

        # get initial parameters
        ypar, yfixed = self.H.ypar_init, self.H.yfixed_inits

        ypar = self.optimize(ypar, yfixed, self._calculate_wavenumbers, niter,
                             derivative, tol, lapack_driver, step_size)

        return ypar

    def optimize_intensity_by_svd(self, niter=1, derivative='n', tol=0.1,
                                  lapack_driver='gelsd', step_size=1.0e-4):

        # get initial parameters
        ypar, yfixed = self.H.ypar_init, self.H.yfixed_init

        # add initial dmf parameters to ypar
        if self.S.dmfs_init is not None:
            for (n, k) in self.S.dmfs_init:
                ypar = np.hstack((ypar, self.S.dmfs_init[(n, k)][:, 1]))
                yfixed = np.hstack((yfixed, self.S.dmfs_init[(n, k)][:, 2]))

        ypar = self.optimize(ypar, yfixed, self._calculate_intensity, niter,
                             derivative, tol, lapack_driver, step_size)

        return ypar

    def get_singular_values(self):
        """Get the singular values generated during the final iteration of the fit

        Returns:
            array: The array of singular values
        """

        return self.sing_values

    def _calculate_numerical_derivatives(self, ypar, yfixed, calculate_data,
                                         ycal_init, step_size):
        """Compute the derivatives matrix by the parameters

        Args:
            ypar (array): fitting parameters
            yfixed (array): fixed/free parameter
            ycal_init (array): the initial calculated energies
            step_size (float): determines the changes in the parameters

        Returns:
            array: the matrix with the derivatives
        """

        dydp = np.zeros(
            shape=(ycal_init.shape[0], ypar.shape[0]), dtype=np.float64)

        # change all free parameters
        delta_par = np.abs(ypar) * step_size
        delta_par[yfixed == 0] = 0.0

        for prm in range(0, ypar.shape[0]):
            if not yfixed[prm]:
                continue
            changed_pars = np.copy(ypar)
            changed_pars[prm] += delta_par[prm]
            ycal, *_ = calculate_data(changed_pars)
            dydp[:, prm] = (ycal - ycal_init) / delta_par[prm]

        return dydp

    def _calculate_analytical_derivatives(self, ypar, yfixed, qnums):

        msize = self.H.nch*self.H.ngrid
        evecs_found = np.zeros((msize, qnums.shape[0]))
        dydp = np.zeros((evecs_found.shape[1], ypar.shape[0]))

        for i in range(0, qnums.shape[0]):
            ename = \
                f'evector_J{int(qnums[i, 1])}_' + \
                f'p{int(qnums[i, 2])}_i{int(qnums[i, 3])}.npy'
            evectors = np.load(_join('eigenvectors', ename))
            vi = np.where(evectors[0, :].astype(int) == int(qnums[i, 0]))[0]
            sti = np.where(evectors[1, :].astype(int) == int(qnums[i, 4]))[0]
            evecs_found[:msize, i] = evectors[2:, np.intersect1d(vi, sti)[0]]

        ejqnums, fjqnums = qnums[qnums[:, 2] == 1], qnums[qnums[:, 2] == 0]
        fjqind = (fjqnums[:, 1] - np.min(fjqnums[:, 1])).astype(int)
        ejqind = (ejqnums[:, 1] - np.min(ejqnums[:, 1])).astype(int)
        ejqind += self.H.jqnumbers.shape[0]
        jqind = np.hstack((fjqind, ejqind))

        for prm in range(0, ypar.shape[0]):
            if not yfixed[prm]:
                continue
            sk = np.zeros((msize, msize))
            for iprm in self.H.block_index[prm]:
                r1, r2, c1, c2, count, cp = iprm
                sk_grid = self.H.sk_grid[r1:r2, count]
                if cp < 0:  # potential
                    np.fill_diagonal(sk[r1:r2, c1:c2], sk_grid)
                    dydp[:, prm] = np.diag(evecs_found.T @ sk @ evecs_found)
                else:  # coupling
                    sk_coef = self.H.interact.cp_map[cp][r1:r2, :]
                    memoize = set()
                    for countj, j in enumerate(jqind):
                        if j not in memoize:
                            elem = self._get_dydp_elem(iprm[:4], sk, sk_grid,
                                                      sk_coef, j, evecs_found)
                        dydp[countj, prm] = elem[countj]
                        memoize.add(j)
        return dydp

    # @staticmethod
    def _get_dydp_elem(self, indices, sk, sk_grid, sk_coef, j, evecs_found):
        r1, r2, c1, c2 = indices
        np.fill_diagonal(sk[r1:r2, c1:c2], sk_grid)
        sk[r1:r2, c1:c2] = sk[r1:r2, c1:c2] * sk_coef[:, j]
        if r1 != c1 or r2 != c2:
            np.fill_diagonal(sk[c1:c2, r1:r2], sk_grid)
            sk[c1:c2, r1:r2] = sk[c1:c2, r1:r2] * sk_coef[:, j]
        elem = np.diag(evecs_found.T @ sk @ evecs_found)

        return elem

    @staticmethod
    def _find_corrections(data_matrix, rhs_vector, tol, lapack_driver):
        """Find corrections to a set of parameters by solving
           a linear system of equations with SVD algorithm

        Args:
            data_matrix (array): The derivative matrix
            rhs_vector (array): RHS vector
            tol (float): Tolerance value
            lapack_driver (str): The lapack routine used for
                nonliner least-squares fiitting

        Returns:
            tuple: the corrections to the parameters, the rank of the
                matrix and an array containing the singular values

        Remarks:
            1. There are two similar procedures for SVD in numpy and
            scipy libraries. In their older versions the definitions of
            the default tol value was different.
            2. Tol controls which and how many linear combinations of
            parameters will be ignored because of the singularity of
            the matrix.
        """

        x, _, rank, s = sp.linalg.lstsq(data_matrix, rhs_vector, tol, lapack_driver)
        # x, _, rank, s = np.linalg.lstsq(A, b, rcond=tol)

        return x, rank, s

    def print_info(self, itr, tol, stats):
        """Print a short information about the progress of the fit on
        each iteration

        Args:
            itr (int): the number of the iteration
            tol (_type_): the tolerance value
            stats (_type_): usefull statistics

        Returns:
            str: The string which to be outputed.
        """

        stats_str = _utils.print_stats(stats)
        out_str = (f'{14*"- "}Iteration {itr} '
                   f'{14*"- "}\nTOL = {tol:.1e}\n')
        out_str += stats_str

        return out_str

    def _format_params_table(self, data):

        out = ''
        init_line = f'  {107*"-"}\n'
        out += init_line
        width = 15
        prec = 6

        title = (f" | {'id':^{width}} | {'Initial value':^{width}} |"
                 f"{'Correction':^{width}} | {'Final value':^{width}}"
                 f" | {'% change':^{width}} | {'Fixed':^{width}} |\n")

        out += title
        out += init_line
        for row in data:
            line = (f' | {row[0]:^{width}.0f} | {row[1]:>{width}.{prec}f}'
                    f' | {row[2]:>{width}.{prec}f} | {row[3]:>{width}.{prec}f}'
                    f' | {row[4]:>{width}.{prec}f} | {row[5]:^{width}.0f} |\n')
            out += line
        out += init_line

        return out

    def _write_fit_progress(self, rank, sing_values, tol, ypar, dypar, yfixed):
        """Write formatted detailed information about the progress of the fit

        Args:
            rank (float): The rank of the SVD matrix
            sing_values (array): An array with singular values
            tol (float): Tolerance value
            ypar (array): Fitting parameters
            dypar (array): Corrections to the parameters
            yfixed (array): Whether a parameter is fixed or free

        Remarks:
            Singular values are sorted from the largest to the smallest one.
            Singular values smaller than tol times the largest singular value
            are set to 0.

            The conditon number is the ratio of the largest to the smallest
            singular value. The larger the condition number the closer to
            singular is the matrix.
        """

        self.progress_str.write(
            f'\nTOL * max singular value = {tol*sing_values[0]:.3f}\n'
            f'\nEffective rank = {rank}'
            f'\nThe first {rank} significant singular values:\n')
        self.progress_str.write(' '.join(
            ['{:.1e}'.format(i) for i in sing_values[:rank]]))

        with np.errstate(divide='ignore', invalid='ignore'):
            cond_number = abs(sing_values[0] / sing_values[-1])

        self.progress_str.write(f'\nCondition number = {cond_number:.1e}\n')

        # TODO: check if there are parameters equal to zero
        pcount = np.arange(1, ypar.shape[0]+1, 1, dtype=np.int64)
        xypar = ypar + dypar
        change = np.absolute(dypar / ypar) * 100

        iypar = ypar  # / self.H.yunits
        corr = dypar  # / self.H.yunits
        fypar = xypar  # / self.H.yunits
        data_out = np.column_stack((pcount, iypar, corr, fypar, change, yfixed))

        self.progress_str.write(self._format_params_table(data_out))

    def get_covariance_matrix(self):
        """Return the covariance matrix
        """
        raise NotImplementedError()
        # find covariance matrix
        # U, s, V = sp.linalg.svd(A, full_matrices=False)
        # si = np.zeros(s.shape[0])
        # print(s)
        # for j in range(s.shape[0]):
        #     if s[j] < max(s)*tol:
        #         si[j] = 0
        #     else:
        #         si[j] = 1.0/(s[j]*s[j])
        # # V = Vt.T
        # cvm = np.zeros((V.shape[0], V.shape[0]))
        # for i in range(V.shape[0]):
        #     for j in range(i+1):  # range(V.shape[0]):
        #         summ = 0.0
        #         for k in range(V.shape[0]):
        #             summ += (V[i, k] * V[j, k]) * si[k]
        #         cvm[j][i] = cvm[i][j] = summ
        # np.savetxt('cov.dat', cvm, fmt='%12.4e')

        # crm = np.zeros((V.shape[0], V.shape[0]))
        # for i in range(s.shape[0]):
        #     for j in range(s.shape[0]):
        #         crm[i, j] = \
        #           cvm[i, j] / (math.sqrt(cvm[i, i]) * math.sqrt(cvm[j, j]))

        # np.savetxt('crm.dat', crm, fmt='%12.4e')

    def optimize_levels_by_svd_regular(self, niter=1, derivative='n',
                                       tol=0.1, lapack_driver='gelsd',
                                       step_size=1.0e-4, regular=0.0):
        pass

    def _minuit_least_squares(self, ypar):
        """The procedure which Minuit calls on each iteration

        Args:
            ypar (numpy ndarray): The fitting parameters

        Returns:
            float: Chi Square value

        Remarks:
            A single argument is required for this function
        """
        ycal, yobs, ydel, yvar, stats = self._calculate_levels(ypar)

        return np.sum(ydel ** 2 / yvar)

    def _set_constraints(self):
        pass

    def optimize_levels_by_minuit(self, niter=0, step_size=1.0e-4, uncert=False):
        """Minuit algorithm for nonliner least squares fitting

        Args:
            niter (int, optional): The number of iterations. Defaults to 0.
            step_size (float, optional): the initial change in the parameters.
                Defaults to 1.0e-4.
            uncert (bool, optional): whether to compute the covariance and the
                correlations matrices and the uncertanties in the parameters.
                Defaults to False.
        """

        print_level = int(self.progress) + 1
        ypar, yfixed = self.H.params, self.H.fixed

        # TODO: # check if this is correct
        errors = np.abs(ypar) * step_size

        m = Minuit.from_array_func(
            self._minuit_least_squares,
            ypar,
            errors,
            fix=1-yfixed,
            limit=None,
            errordef=Minuit.LEAST_SQUARES,
            print_level=print_level
        )

        print(f'Running Minuit Fit...\nInitial parameters:\n{m.params}\n')

        m.migrad(ncall=niter)

        self.progress_str = io.StringIO()

        self.progress_str.write(
            f'Minuit Fit completed.\n'
            f'The initial parameters were:\n{m.init_params}\n'
            f'The final prameters are:\n{m.params}\n'
            f'The reduced chi2 = {m.fval / (len(ypar) - m.nfit)}'
        )
        # reduced chi2 - should be around 1
        # print(m.fmin)
        # print(m.migrad_ok())
        # print(m.fval)
        # print(m.matrix_accurate())

        # TODO: check how to get the uncertenties
        if uncert:
            np.set_printoptions(linewidth=150)
            np.set_printoptions(precision=3)
            self.progress_str.write(
                f'\nCorrelation matrix\n{m.np_matrix(correlation=True)}\n'
                f'\nCovariance matrix\n{m.np_covariance()}\n'
                f'\nParameter uncertainties\n{m.hesse()}\n'
            )

        print(self.progress_str.getvalue())

    def run_levmar(self, niter=1):

        """Levenberg-Marquardt algorithm for nonliner least squares fitting

        Args:
            niter (int, optional): Number of iterations. Defaults to 1.
        """

        raise NotImplementedError('To be implemented!')
