import os
import io
import math
import numpy as np
import scipy as sp
from .molecule_data import Channel, Coupling
from .utils import Utils, C_hartree
# from numba import njit, guvectorize, int64, float64

try:
    from iminuit import Minuit
except ModuleNotFoundError:
    pass
    # print("'iminuit' module is not installed!\n")

__all__ = ['Fitting']


class Fitting:

    def __init__(self, ml, progress=False):

        self.ml = ml
        self.progress = progress

        try:
            Utils.createBackup(Coupling.cpl_file)
        except AttributeError:
            pass

    def _get_initial_parameters(self):

        if len(self.ml.couplings) > 0:
            ypar = np.hstack((Channel.ppar, Coupling.cpar))
            yfixed = np.hstack((Channel.pfixed, Coupling.cfixed))
            preg = np.zeros(Channel.ppar.shape[0])
            plambda = np.zeros(Channel.ppar.shape[0])
            yreg = np.hstack((preg, Coupling.cregular))
            ylam = np.hstack((plambda, Coupling.clambda))
        else:
            ypar = np.copy(Channel.ppar)
            yfixed = np.copy(Channel.pfixed)
            yreg = np.zeros(Channel.ppar.shape[0])
            ylam = np.zeros(Channel.ppar.shape[0])

        return ypar, yfixed, yreg, ylam

    def _get_eigenvalues(self, ypar, is_weighted=False):

        self.ml.store_predicted = False
        self.ml.is_weighted = is_weighted
        stats = self.ml.calculate_eigenvalues(ypar=ypar)

        # Uncomment this if intermediate results are needed for debugging
        print(f'Chi Square = {stats[0]:<18.8f} |RMS = {stats[1]:<15.8f}cm-1')

        return self.ml.out_data, stats

    def _minuit_least_squares(self, ypar):

        """The procedure which Minuit calls on each iteration

        Args:
            ypar (numpy ndarray): The fitting parameters

        Returns:
            float: Chi Square value

        Remarks:
            1. Minuit requires this procedure to accept a single argument
        """

        out_data, _ = self._get_eigenvalues(ypar)
        ydel = out_data[:, 9] / C_hartree
        yvar = out_data[:, 10] / C_hartree

        Channel.edit_channel_parameters(ypar, self.ml.channels)

        if len(self.ml.couplings) > 0:
            Coupling.edit_coupling_parameters(ypar, self.ml.couplings)

        return np.sum(ydel ** 2 / yvar)

    def set_constraints(self):
        pass

    def run_minuit(self, niter=0, step_size=1.0e-4,
                   set_limits=None, uncert=False):
        """Minuit algorithm for nonliner least squares fitting

        Args:
            niter (int, optional): The number of iterations. Defaults to 0.
            step_size (float, optional): the initial change in the parameters.
            Defaults to 1.0e-4.
            set_limits (array, optional): limit the possible values of
            the parameters. Defaults to None.
            uncert (bool, optional): whether to compute the covariance and the
            correlations matrices and the uncertanties in the parameters.
            Defaults to False.
        """

        print_level = int(self.progress) + 1

        ypar, yfixed, *_ = self._get_initial_parameters()

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

    def run_svd(self, niter=1, deriv='n', tol=0.1, lapack_driver='gelsd',
                step_size=1.0e-4, is_weighted=False, regular=0.0):

        """SVD algorithm for nonliner least squares fitting.Solves Ax = b using
        Singular Value Decomposition (SVD). Starting from some trial set of
        parameters on each iteration the linear system Ax=b for the corrections
        x to the current set of parameters x_cur will be solved and a new set
        x_new = x_cur + x will be obtained. Matrix A contains the derivatives
        of the energies with respect to the parameters and b is the rhs vector
        E_obs - E_cal(x_cur).

        Args:
            niter (int, optional): number of iterations. Defaults to 1.
            deriv (str, optional): Type of derivative. Defaults to 'n'.
            tol (float, optional): the tolerance - determines which linear
            combinations of parameters to be discarded by SVD. Defaults to 0.1.
            lapack_driver (str, optional): lapack routine. Defaults to 'gelsd'.
            step_size (float, optional): the value used to change the
            parameters in the calculation of derivatives. Defaults to 1.0e-4.
            is_weighted (bool, optional): weighted fit. Defaults to False.
            restart (bool, optional): Not yet implemented. Defaults to False.
        """

        # (0) get the initial parameters
        ypar, yfixed, yreg, ylam = self._get_initial_parameters()
        is_failed = False

        for it in range(1, niter + 1):

            self.progress_str = io.StringIO()
            self.progress_str.write(
                f'\n{18*"- "}Iteration {it} {18*"- "} \n'
                f'\nTOL{"":>23}= {tol:.1e}'
            )

            # (1) get initially calculated energies on every new iteration
            out_data, stats = self._get_eigenvalues(
                ypar, is_weighted=is_weighted
            )

            ycal = out_data[:, 7].astype(np.float64) / C_hartree
            ydel = -out_data[:, 9].astype(np.float64) / C_hartree
            yvar = out_data[:, 10].astype(np.float64) / C_hartree

            chi2_best, rms_best = stats[0], stats[1]

            if regular != 0.0:
                ydiff = yreg - ypar
                cmult = np.divide(
                    ydiff, ylam, out=np.zeros_like(ydiff), where=ylam != 0
                )
                c = cmult * regular
                chi2r = np.sum(np.square(c))
                chi2_best = chi2_best + chi2r
                rms_best = rms_best + math.sqrt(chi2r)

            self.progress_str.write(
                f'\nChi square initial{"":>8}= {chi2_best:.8f}'
                f'\nRMS initial{"":>15}= {rms_best:.8f} cm-1\n'
            )

            # (2) generate the matrix A - the matrix containing the
            # derivatives of the energies with respect to the parameters

            if deriv == 'n':
                dydp = self._generate_numerical_derivatives(
                    ypar, yfixed, ycal, is_weighted=is_weighted,
                    step_size=step_size
                )
                # np.savetxt('dydp_numerical.dat', dydp, fmt='%.6e')
            else:
                # the quantum numbers are aranged as follows:
                # v, J, par, iso, state
                qnums = np.column_stack((
                    out_data[:, 0], out_data[:, 1], out_data[:, 5],
                    (out_data[:, 6] % self.ml.nisotopes)+1, out_data[:, -1]
                ))
                dydp = self._generate_analytical_derivatives(
                    ypar, yfixed, qnums, is_weighted=is_weighted
                )
                # np.savetxt('dydp_analytical.dat', dydp, fmt='%.6e')

            A = (dydp.T / yvar).T

            # (3) calculate the rhs vector b
            b = ydel / yvar

            # extend the system Ax=b with new one: Bx=c
            if regular != 0.0:
                # B = (np.eye(ypar.shape[0], ypar.shape[0]) / ylam) * regular
                idm = np.eye(ypar.shape[0], ypar.shape[0])
                bmult = np.divide(
                    idm, ylam, out=np.zeros_like(idm), where=ylam != 0
                )
                B = bmult * regular
                A = np.vstack((A, B))
                ydiff = yreg - ypar
                cmult = np.divide(
                    ydiff, ylam, out=np.zeros_like(ydiff), where=ylam != 0
                )
                c = cmult * regular
                b = np.hstack((b, c))

            # the proposed corrections to the parameters is the solution
            x, rank, sv = Fitting._find_corrections(A, b, tol, lapack_driver)

            # (4) calculate new energies using the proposed corrections
            _, stats_try = self._get_eigenvalues(
                ypar+x, is_weighted=is_weighted
            )

            chi2_try, rms_try = stats_try[0], stats_try[1]

            if regular != 0.0:
                chi2r = np.sum(np.square((c - np.dot(B, x))))
                chi2_try = chi2_try + chi2r
                rms_try = rms_try + math.sqrt(chi2r)

            # eps = 1.0e-5

            # (5) Try to change the corrections
            is_failed = False

            if chi2_try < chi2_best:  # abs(chi2_try - chi2_best) < eps:
                chi2_best = chi2_try
                rms_best = rms_try
                ypar = ypar + x
            else:  # chi2try >= chi2best
                step = 10.0
                ypar_try = np.copy(ypar)
                x_try = np.copy(x)
                x_try = x_try / step

                # try to improve the chi_square value
                while True:
                    # x_try = x_try / step
                    ypar_try = ypar_try + x_try

                    _, stats_new = self._get_eigenvalues(
                        ypar_try, is_weighted=is_weighted
                    )

                    chi2_new, rms_new = stats_new[0], stats_new[1]

                    if regular != 0.0:
                        chi2r = np.sum(
                            np.square(((yreg-ypar_try)) - np.dot(B, x_try))
                        )
                        chi2_new = chi2_new + chi2r
                        rms_new = rms_new + math.sqrt(chi2r)

                    # if abs(chi2_new - chi2_best) < eps:
                    if chi2_new < chi2_best:
                        chi2_best = chi2_new
                        rms_best = rms_new
                        ypar = ypar_try

                    else:  # chi2new >= chi2best:
                        is_failed = True

                        self.progress_str.write(
                            f'\nCorrections changed =>\nChi square = '
                            f'{chi2_best:.8f}\n'
                            f'RMS = {rms_best:.8f}\n'
                        )

                        break

            # if it isn't got improved then change tol
            if is_failed:
                tol /= 5.0
                self.progress_str.write(f'\nTOL CHANGED {"":>15}= {tol:.1e}')

            # (6) Save the final parameters
            if self.progress:
                self._save_svd_parameters(rank, sv, tol, ypar, x, yfixed, it)

            Channel.edit_channel_parameters(ypar, self.ml.channels)

            if len(self.ml.couplings) > 0:
                Coupling.edit_coupling_parameters(ypar, self.ml.couplings)

            self.progress_str.write(
                f'\nChi square final/best{"":>6}= {chi2_best:.8f} '
                f'for iteration {it}'
                f'\nRMS final{"":>18}= {rms_best:.8f} cm-1'
            )

            print(self.progress_str.getvalue())

    def _generate_numerical_derivatives(self, ypar, yfixed, ycal_init,
                                        is_weighted, step_size):
        """Compute the derivatives matrix by the parameters

        Args:
            ypar (array): fitting parameters
            yfixed (array): fixed/free parameter
            ycal_init (array): the initial calculated energies
            is_weighted (bool): whether to apply weighted fit
            step_size (float): determines the changes in the parameters

        Returns:
            array: the matrix with the derivatives
        """

        dydp = np.zeros(
            shape=(ycal_init.shape[0], ypar.shape[0]), dtype=np.float64
        )
        # change all free parameters
        dp = np.abs(ypar) * step_size
        dp[yfixed == 0] = 0.0

        for prm in range(0, ypar.shape[0]):
            if yfixed[prm]:
                dpar = np.copy(ypar)
                dpar[prm] += dp[prm]
                out_data, _ = \
                    self._get_eigenvalues(dpar, is_weighted=is_weighted)
                ycal = out_data[:, 7] / C_hartree
                dydp[:, prm] = (ycal - ycal_init) / dp[prm]

        return dydp

    def _generate_analytical_derivatives(self, ypar, yfixed,
                                         qnums, is_weighted):

        if not self.ml.sderiv:
            raise SystemExit(
                'Error: Spline derivatives should be calculated. '
                'Change model to cspline.'
            )

        size = self.ml.nch*self.ml.ngrid
        evec_selected = np.zeros((size, qnums.shape[0]))

        for i in range(0, qnums.shape[0]):
            evec_name = \
                f'evector_J{int(qnums[i, 1])}_' + \
                f'p{int(qnums[i, 2])}_i{int(qnums[i, 3])}.npy'

            evectors = np.load(os.path.join('eigenvectors', evec_name))
            vibnums = evectors[0, :].astype(np.int64)
            states = evectors[1, :].astype(np.int64)
            vi = np.where(vibnums == int(qnums[i, 0]))[0]
            sti = np.where(states == int(qnums[i, 4]))[0]
            evector = evectors[2:, np.intersect1d(vi, sti)[0]]
            evec_selected[0:size, i] = evector

        dydp = np.zeros((evec_selected.shape[1], ypar.shape[0]))

        for prm in range(0, ypar.shape[0]):
            if yfixed[prm]:
                sk = np.zeros((size, size))
                for ind in self.ml.block_index[prm]:
                    r1, r2, c1, c2, count = ind
                    if r1 == c1 and r2 == c2:
                        # not working for non-diagonal block
                        sk = sk + np.diag(self.ml.sk_grid[:, count])
                    else:
                        # i don't know why this is not working in both cases
                        np.fill_diagonal(
                            sk[r1:r2, c1:c2], self.ml.sk_grid[:, prm]
                        )
                        np.fill_diagonal(
                            sk[c1:c2, r1:r2], self.ml.sk_grid[:, prm]
                        )
                        # sk = sk + sk.T  # <-- not working for diagonal block
                enr_vec = np.diag(evec_selected.T @ sk @ evec_selected)
                dydp[:, prm] = enr_vec

        return dydp

    @staticmethod
    def _find_corrections(A, b, tol, lapack_driver):
        """Find the corrections to the parameters by solving
        the system of linear equations with SVD technique

        Args:
            A (array): matrix with the derivatives
            b (array): rhs vector
            tol (float): the tolerance value
            lapack_driver (str): the lapack routine

        Returns:
            tuple: the corrections, the rank of the
            matrix and the matrix with the singular values

        Remarks:
            1. There are two similar procedures for SVD from numpy
                and scipy. In their older versions the definitions of
                the default tol value is different
            2. Tol controls which and how many linear combinations of the
                parameters will be ignored since the matrix is singular
        """

        x, _, rank, s = sp.linalg.lstsq(A, b, tol, lapack_driver)
        # x, _, rank, s = np.linalg.lstsq(A, b, rcond=tol)

        return x, rank, s

    def _save_svd_parameters(self, rank, sv, tol, ypar, x, yfixed, it):

        # singular values are sorted from the largest to the smallest one
        # singular values smaller than tol*largest_singular_value are set to 0

        self.progress_str.write(
            f'\nTOL * max singular value{"":>2}= {tol*sv[0]:.3f}\n'
            f'\nSVD matrix rank = the number significant singular values '
            f'= {rank}'
            f'\nThe first {rank} significant singular values:\n'
        )

        self.progress_str.write(
            ' '.join(['{:.1e}'.format(i) for i in sv[0:rank]])
        )

        # the conditon number is the ratio of the largest to the
        # smallest singular value; the larger the condition number
        # the closer to singular is the matrix
        # if cond. number = inf => the matrix is singular

        with np.errstate(divide='ignore', invalid='ignore'):
            cond_number = abs(sv[0] / sv[-1])

        self.progress_str.write(
            f'\nCondition number{"":>10}= {cond_number:.1e}'
        )
        self.progress_str.write(
            f'\n\n{21*"- "}Parameters after iteration {it} {20*"- "}'
        )

        header = (
            f'\n{111*"-"}\n|{"":>9}|{"":>5}Initial value{"":>4}|'
            f'{"":>5}Correction{"":>7}|{"":>5}Final value{"":>6}|'
            f'{"":>5}Percent change{"":>3}|{"":>1}Fixed{"":>1}|\n{111*"-"}\n'
        )
        self.progress_str.write(header)

        # TODO: what if some parameter has 0 value?
        nn = np.arange(1, ypar.shape[0]+1, 1, dtype=np.int64)
        xypar = ypar + x
        change = np.absolute(x/ypar) * 100

        params_str = ''
        for i in range(0, nn.shape[0]):
            params_str = (
                f'|{"":>1}{i:4d}{"":>4}|{"":>3}{ypar[i]:>16.12f}{"":>3}|'
                f'{"":>3}{x[i]:>16.12f}{"":>3}|{"":>3}'
                f'{xypar[i]:>16.12f}{"":>3}|'
                f'{"":>3}{change[i]:>16.12f}{"":>3}|'
                f'{"":>1}{yfixed[i]:>3d}{"":>3}|\n'
            )

            self.progress_str.write(params_str)

    def run_levmar(self, niter=1):

        """Levenberg-Marquardt algorithm for nonliner least squares fitting

        Args:
            niter (int, optional): Number of iterations. Defaults to 1.
        """

        pass
