import io
import math
from os.path import join as _join
from math import sqrt as _sqrt
from logging import error as _error, warning as _warning
import numpy as np
import scipy as sp
from molecule_data import Channel, Coupling
from utils import Utils, C_hartree
# from numba import njit, guvectorize, int64, float64

try:
    from iminuit import Minuit
except ModuleNotFoundError:
    pass
    # print("'iminuit' module is not installed!\n")

__all__ = ['Fitting']


class Fitting:

    def __init__(self, ml, progress_full=False):

        self.ml = ml
        self.progress = progress_full

        try:
            Utils.createBackup(Coupling.cpl_file)
        except AttributeError:
            pass

    def _get_initial_parameters(self):

        if len(self.ml.couplings) > 0:
            ypar = np.hstack((Channel.ppar, Coupling.cpar))
            yunits = np.hstack((Channel.punits, Coupling.cunits))
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
        # print(yunits)
        # print(ypar.shape, yunits.shape)
        return ypar, yfixed, yreg, ylam

    def _get_eigenvalues(self, ypar, is_weighted=False):

        self.ml.store_predicted = False
        self.ml.is_weighted = is_weighted
        stats = self.ml.calculate_eigenvalues(ypar=ypar)

        # Uncomment this if intermediate results are needed for debugging
        # print(f'Chi Square = {stats[0]:<18.8f} |RMS = {stats[1]:<15.8f}cm-1')

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

        deriv = deriv.lower()
        if deriv == 'a':
            self.check_derivative()

        # (0) get the initial parameters
        ypar, yfixed, yreg, ylam = self._get_initial_parameters()
        eps = 1.0e-3
        chi2_best, rms_best = 0.0, 0.0

        print(f'Total number of parameters = {ypar.shape[0]}')
        print(
            f'Total number of optimized parameters = '
            f'{yfixed[yfixed == 1].shape[0]}\n'
        )

        for it in range(1, niter + 1):

            # (1) get the updated energies from the previous iteration
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
                rms_best = rms_best + _sqrt(chi2r)

            info_init = (
                f'{14*"- "}Iteration {it} {14*"- "}\nTOL = {tol:.1e}\n'
                f'Chi square = {chi2_best:.6f} | RMS = {rms_best:.6f} cm-1\n'
            )
            print(info_init)
            self.progress_str = io.StringIO()
            self.progress_str.write(info_init)

            # (2) generate the matrix containing the derivatives of the
            # energies with respect to the parameters

            if deriv == 'n':
                dydp = self._generate_numerical_derivatives(
                    ypar, yfixed, ycal, is_weighted=is_weighted,
                    step_size=step_size
                )
                # np.savetxt('dydp_numerical.dat', dydp, fmt='%.6e')
            else:
                # v, J, par, iso, state
                qnums = np.column_stack((
                    out_data[:, 0], out_data[:, 1], out_data[:, 5],
                    (out_data[:, 6] % self.ml.nisotopes)+1, out_data[:, -1]
                ))
                dydp = self._generate_analytical_derivatives(
                    ypar, yfixed, qnums
                )
                # np.savetxt('dydp_analytical.dat', dydp, fmt='%.6e')

            # (3) calculate the data matrix A from derivatives matrix
            # then calculate the rhs vector b
            A = (dydp.T / yvar).T
            b = ydel / yvar

            # extend the system Ax=b with new one: Bx=c
            if regular != 0.0:
                # B = (np.eye(ypar.shape[0], ypar.shape[0]) / ylam) * regular
                identm = np.eye(ypar.shape[0], ypar.shape[0])
                bmult = np.divide(
                    identm, ylam, out=np.zeros_like(identm), where=ylam != 0
                )
                B = bmult * regular
                A = np.vstack((A, B))
                ydiff = yreg - ypar
                cmult = np.divide(
                    ydiff, ylam, out=np.zeros_like(ydiff), where=ylam != 0
                )
                c = cmult * regular
                b = np.hstack((b, c))

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
            #         crm[i, j] = cvm[i, j] / (math.sqrt(cvm[i, i]) * math.sqrt(cvm[j, j]))

            # np.savetxt('crm.dat', crm, fmt='%12.4e')
            # the propossed corrections to the parameters is the solution
            x, rank, sv = Fitting._find_corrections(A, b, tol, lapack_driver)

            # (4) calculate new energies using the proposed corrections
            _, stats_try = self._get_eigenvalues(
                ypar+x, is_weighted=is_weighted
            )
            chi2_try, rms_try = stats_try[0], stats_try[1]

            if regular != 0.0:
                chi2r = np.sum(np.square((c - np.dot(B, x))))
                chi2_try = chi2_try + chi2r
                rms_try = rms_try + _sqrt(chi2r)

            # (5) Try to change the corrections
            if chi2_try + eps < chi2_best:
                chi2_best, rms_best = chi2_try, rms_try
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

                    _, stats_new = self._get_eigenvalues(
                        ypar_try, is_weighted=is_weighted
                    )
                    chi2_new, rms_new = stats_new[0], stats_new[1]

                    if regular != 0.0:
                        chi2r = np.sum(
                            np.square(((yreg-ypar_try)) - np.dot(B, x_try))
                        )
                        chi2_new = chi2_new + chi2r
                        rms_new = rms_new + _sqrt(chi2r)

                    if chi2_new + eps < chi2_best:
                        chi2_best, rms_best = chi2_new, rms_new
                        ypar = ypar_try
                        info_corr = (
                            f'Smaller step => Chi square = '
                            f'{chi2_best:.6f} | RMS = {rms_best:.6f}'
                        )
                        print(info_corr)
                        self.progress_str.write(info_corr)
                        is_failed = False
                    else:  # chi2new >= chi2best:
                        break

            # on odd iteration try to change tol if nothing has been improved
            if it > 1 and it % 2 == 1 and is_failed:
                tol /= 5.0
                info_tol = f'\nTol changed = {tol:.1e}'
                print(info_tol)
                self.progress_str.write(info_tol)

            # (6) Save the final parameters
            if self.progress:
                self._save_svd_parameters(rank, sv, tol, ypar, x, yfixed, it)

            Channel.edit_channel_parameters(ypar, self.ml.channels)

            if len(self.ml.couplings) > 0:
                Coupling.edit_coupling_parameters(ypar, self.ml.couplings)

            if self.progress:
                print(self.progress_str.getvalue())

        self.singular_values = sv

        print(
            f'\nChi square final = {chi2_best:.6f} | '
            f'RMS final = {rms_best:.6f} cm-1'
        )

    def get_singular_values(self):

        return self.singular_values

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
            if not yfixed[prm]:
                continue
            dpar = np.copy(ypar)
            dpar[prm] += dp[prm]
            out_data, _ = self._get_eigenvalues(dpar, is_weighted=is_weighted)
            ycal = out_data[:, 7] / C_hartree
            dydp[:, prm] = (ycal - ycal_init) / dp[prm]

        return dydp

    def _generate_analytical_derivatives(self, ypar, yfixed, qnums):

        msize = self.ml.nch*self.ml.ngrid
        evecs_found = np.zeros((msize, qnums.shape[0]))
        dydp = np.zeros((evecs_found.shape[1], ypar.shape[0]))

        for i in range(0, qnums.shape[0]):
            ename = \
                f'evector_J{int(qnums[i, 1])}_' + \
                f'p{int(qnums[i, 2])}_i{int(qnums[i, 3])}.npy'
            evectors = np.load(_join('eigenvectors', ename))
            vi = np.where(evectors[0, :].astype(int) == int(qnums[i, 0]))[0]
            sti = np.where(evectors[1, :].astype(int) == int(qnums[i, 4]))[0]
            # print(qnums[i])
            # print(vi)
            # print(sti)
            # print(evectors[2:, np.intersect1d(vi, sti)])
            evecs_found[:msize, i] = evectors[2:, np.intersect1d(vi, sti)[0]]

        ejqnums, fjqnums = qnums[qnums[:, 2] == 1], qnums[qnums[:, 2] == 0]
        fjqind = (fjqnums[:, 1] - np.min(fjqnums[:, 1])).astype(int)
        ejqind = (ejqnums[:, 1] - np.min(ejqnums[:, 1])).astype(int)
        ejqind += self.ml.jqnumbers.shape[0]
        jqind = np.hstack((fjqind, ejqind))

        for prm in range(0, ypar.shape[0]):
            if not yfixed[prm]:
                continue
            sk = np.zeros((msize, msize))
            for iprm in self.ml.block_index[prm]:
                r1, r2, c1, c2, count, cp = iprm
                sk_grid = self.ml.sk_grid[r1:r2, count]
                if cp < 0:  # potential
                    np.fill_diagonal(sk[r1:r2, c1:c2], sk_grid)
                    dydp[:, prm] = np.diag(evecs_found.T @ sk @ evecs_found)
                else:  # coupling
                    sk_coef = self.ml.interact.cp_map[cp][r1:r2, :]
                    memoize = set()
                    for countj, j in enumerate(jqind):
                        if j not in memoize:
                            elem = self.get_dydp_elem(iprm[:4], sk, sk_grid,
                                                      sk_coef, j, evecs_found)
                        dydp[countj, prm] = elem[countj]
                        memoize.add(j)
        return dydp

    # @staticmethod
    def get_dydp_elem(self, indices, sk, sk_grid, sk_coef, j, evecs_found):
        r1, r2, c1, c2 = indices
        np.fill_diagonal(sk[r1:r2, c1:c2], sk_grid)
        sk[r1:r2, c1:c2] = sk[r1:r2, c1:c2] * sk_coef[:, j]
        if r1 != c1 or r2 != c2:
            np.fill_diagonal(sk[c1:c2, r1:r2], sk_grid)
            sk[c1:c2, r1:r2] = sk[c1:c2, r1:r2] * sk_coef[:, j]
        elem = np.diag(evecs_found.T @ sk @ evecs_found)

        return elem

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

    def check_derivative(self):

        if not self.ml.is_cspline_pec and not self.ml.is_cspline_cpl:
            _error(
                'Change model parameter to cspline: analytical '
                'derivatives can be calculated only for model=cspline'
            )
            raise SystemExit()

        if not self.ml.is_cspline_pec:
            _warning(
                'The parameter model for the channel objects is not set '
                'to cspline.\nThis will cause a failure (KeyError exception) '
                'if you try to optimize the channel parameters.\n'
                'Change model=cspline.\n'
            )

        if not self.ml.is_cspline_cpl:
            _warning(
                'The parameter model for the coupling objects is not set '
                'to cspline.\nThis will cause a failure (KeyError exception) '
                'if you try to optimize the coupling parameters.\n'
                'Change model=cspline.\n'
            )

    def _save_svd_parameters(self, rank, sv, tol, ypar, x, yfixed, it):

        # singular values are sorted from the largest to the smallest one
        # singular values smaller than tol*largest_singular_value are set to 0

        self.progress_str.write(
            f'\nTOL * max singular value{"":>2}= {tol*sv[0]:.3f}\n'
            f'\nEffective rank = {rank}'
            f'\nThe first {rank} significant singular values:\n'
        )

        self.progress_str.write(
            ' '.join(['{:.1e}'.format(i) for i in sv[:rank]])
        )

        # the conditon number is the ratio of the largest to the
        # smallest singular value; the larger the condition number
        # the closer to singular is the matrix

        with np.errstate(divide='ignore', invalid='ignore'):
            cond_number = abs(sv[0] / sv[-1])

        self.progress_str.write(
            f'\nCondition number = {cond_number:.1e}'
        )

        header = (
            f'\n{111*"-"}\n|{"":>8}|{"":>3}Initial value{"":>5}|'
            f'{"":>4}Correction{"":>6}|{"":>4}Final value{"":>6}|'
            f'{"":>1}Percent change{"":>1}|{"":>1}Fixed{"":>1}|\n{111*"-"}\n'
        )
        self.progress_str.write(header)

        # TODO: what if some parameter has 0 value?
        nn = np.arange(1, ypar.shape[0]+1, 1, dtype=np.int64)
        xypar = ypar + x
        change = np.absolute(x/ypar) * 100

        params_str = ''
        for i in range(1, nn.shape[0]):
            params_str = (
                f'|{"":>1}{i:4d}{"":>3}|{"":>2}{ypar[i]:>16.8f}{"":>3}|'
                f'{"":>2}{x[i]:>16.8f}{"":>2}|{"":>3}'
                f'{xypar[i]:>16.8f}{"":>2}|'
                f'{"":>2}{change[i]:>12.8f}{"":>2}|'
                f'{"":>1}{int(yfixed[i]):>3d}{"":>3}|\n'
            )

            self.progress_str.write(params_str)

    def run_levmar(self, niter=1):

        """Levenberg-Marquardt algorithm for nonliner least squares fitting

        Args:
            niter (int, optional): Number of iterations. Defaults to 1.
        """

        raise NotImplementedError()
