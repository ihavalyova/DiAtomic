import os
import io
import time
import math
import numpy as np
import scipy as sp

from pprint import pprint
from sklearn.metrics import mean_squared_error
from more_itertools import unique_everseen

from utils import Utils
from molecule_data import Channel
from molecule_data import Coupling
from constants import Const

class Fitting:

    def __init__(self, ml, overwrite=True, progress_fit=False):

        self.ml = ml
        self.overwrite = overwrite
        self.progress_fit = progress_fit

        try:
            Utils.createBackup(Coupling.cpl_file)
        except AttributeError:
            pass

        for pname in Channel.pnames:
            Utils.createBackup(pname)

    def run(self, method='svd', niter=1, tol=0.1, deriv='n', is_weighted=False):

        """The fit function called by the user for running the fit

        Args:
            method (str, optional): The fit method.
            The options are 'svd', 'minuit', 'levmar'. Defaults to 'svd'.
            niter (int, optional): Number of fit iterations. Defaults to 1.
            tol (float, optional): Tolerance value used only with SVD method. Defaults to 0.1.
            deriv (str, optional): Numerical or Analytical derivatives. Defaults to 'n'.
            is_weighted (bool, optional): Perform weighed fit by Watson's method. Defaults to False.

        Raises:
            SystemExit: if incorrect value for method parameter is provided
        """

        # if deriv != 'n' or deriv != 'a':
        #     deriv = 'n'

        ypar, yfixed = self.get_initial_parameters()

        if method.lower() == 'minuit':
            try:
                from iminuit import Minuit
            except ImportError as err:
                SystemExit(err)

            fit = MinuitFit(self.ml, self.overwrite, self.progress_fit)
            fit.run_minuit(ypar, yfixed, niter)

        elif method.lower() == 'svd':
            fit = SvdFit(self.ml, self.overwrite, self.progress_fit)
            fit.run_svd(ypar, yfixed, niter, deriv, tol, is_weighted)

        elif method.lower() == 'levmar':
            pass
        else:
            raise SystemExit(f'\nFitting procedure {method} is not supported')

    def get_initial_parameters(self):

        """Combine all parameters in one array which to be used in the fit

        Returns:
            tuple of 2 arrays: values of all parameters and which of them is fixed/free
        """

        #ppar, pfixed = Channel.get_channel_parameters()
        ppar, pfixed = Channel.get_unique_channels_parameters(self.ml.channels)
        cpar, cfixed = np.array([]), np.array([])

        if len(self.ml.couplings) > 0:
            cpar, cfixed = Coupling.get_coupling_parameters(self.ml.couplings)

        ypar = np.concatenate((ppar, cpar), axis=None)
        yfixed = np.concatenate((pfixed, cfixed), axis=None)

        return ypar, yfixed


class Fit:

    def __init__(self, ml, overwrite, progress_fit):

        self.ml = ml
        self.overwrite = overwrite
        self.progress_fit = progress_fit

    # TODO: change the name of this procedure
    def getEigenvalues(self, ypar, is_weighted=False):

        #self.ml.calculateLevels(store_predicted=False, ypar=ypar, is_weighted=is_weighted)
        self.ml.store_predicted = False ## ????
        self.ml.is_weighted = is_weighted ## ????
        self.ml.calculate_eigenvalues(ypar=ypar)

        return self.ml.out_data

class MinuitFit(Fit):
    
    """Class defining the Minuit algorithm for nonliner least squares fitting

    Args:
        Fit (object): The inherited class
    """

    def __init__(self, ml, overwrite, progress_fit):

        super().__init__(ml, overwrite, progress_fit)

        self.count_it = 1

    def minuit_least_squares(self, ypar):

        """The procedure which Minuit will call on each iteration

        Args:
            ypar (numpy ndarray): The fitting parameters

        Returns:
            float: Chi Square value

        Remarks:
            1. Minuit requires this procedure to accept a single argument
        """

        #self.edit_yml_file(ypar)

        out_data = super().getEigenvalues(ypar)

        ycal = out_data[:,8]
        yexp = out_data[:,9]
        ydel = out_data[:,10]
        yvar = out_data[:,11]

        #rms = self.calculate_rms(yexp, ycal)
        #rmsd = self.calculate_dimless_rms(yexp, ycal, yvar, is_weighted=False)

        chi2, rms, rmsd = self.ml.calculate_stats(yexp, ycal, yvar, is_weighted=False)

        if self.progress_fit:
            print(f'Chi2 = {chi2:.8f} for iteration {self.count_it}')
            print(f'RMS = {rms:.8f} cm-1 for iteration {self.count_it}')
            print(f'RMSD = {rmsd:.8f} for iteration {self.count_it}\n')

        self.count_it += 1

        if self.overwrite:
            Channel.edit_channel_parameters(ypar, self.ml.channels)

            if len(self.ml.couplings) > 0:
                Coupling.edit_coupling_parameters(ypar, self.ml.couplings)

        return np.sum(ydel ** 2 / yvar)

    def set_constraints(self):
        pass

    def run_minuit(self, ypar, yfixed, niter):

        """The main procedure which is calling the Minuit procedure

        Args:
            ypar (ndarray): fitting parameters
            yfixed (ndarray): a binary array defining whether a parameter is free
            niter (int): the number of fit iterations
        """

        errors = np.abs(ypar) * 1.0e-4

        m = Minuit.from_array_func(
            self.minuit_least_squares, 
            ypar, 
            errors,                       
            fix=yfixed, 
            limit=None,
            errordef=1, 
            print_level=2
        )

        m.print_param()
        fmin, param = m.migrad(ncall=niter)

        pprint(fmin)
        pprint(param)
        m.print_param()
        pprint(m.migrad_ok())
        pprint(m.fval)
        print(m.matrix_accurate())
        pprint(m.np_matrix())
        pprint(m.np_matrix(correlation=True))


class SvdFit(Fit):

    """Class defining the SVD algorithm for nonliner least squares fitting

    Args:
        Fit (object): The inherited class
    """

    def __init__(self, ml, overwrite, progress_fit):

        super().__init__(ml, overwrite, progress_fit)
        self.progress_str = None

    def run_svd(self, ypar, yfixed, niter, deriv, tol, is_weighted):

        """The main procedure for the svd fit

        Args:
            ypar (ndarray): The fitting parameters
            yfixed (ndarray): Which parameter is fixed/free
            niter (int): the number of fit iterations
            tol (float): The tolerance value
            is_weighted (bool): whether to apply weighted fitting 
        """
    
        for it in range(1, niter + 1):

            self.progress_str = io.StringIO()

            self.progress_str.write(f'\n {15*"="} Iteration {it} {15*"="} \n')

            # 1) Get initially calculated data
            out_data = super().getEigenvalues(ypar, is_weighted=is_weighted)

            ycal = out_data[:,8].astype(np.float64)
            yexp = out_data[:,9].astype(np.float64)
            ydel = -out_data[:,10].astype(np.float64)
            yvar = out_data[:,11].astype(np.float64)

            chisq_best, rms_init, rmsd_init = self.ml.calculate_stats(yexp, ycal, yvar, is_weighted)

            self.progress_str.write(
                f'\nChi square initial = {chisq_best:.8f} for iteration {it}'
                f'\nRMS initial = {rms_init:.8f} cm-1 for iteration {it}'
                f'\nRMSD initial = {rmsd_init:.8f} for iteration {it}\n'
            )

            # 2) Generate design matrix A
            # A is found through calculation of numerical derivatives
            
            if deriv == 'n':
                dydp = self.generate_design_matrix_numerical(
                    ypar, yfixed, ycal, is_weighted=is_weighted
                )

                #np.savetxt('true_dydp.dat', dydp, fmt='%.6e')
            else:
                dydp = self.generate_design_matrix_analytical(
                    ypar, yfixed, ycal, is_weighted=is_weighted
                )

            # 3) Solve Ax = b using Singular Value Decomposition

            # find rhs vector
            b = ydel / yvar

            # the proposed corrections to the parameters is the solution
            x = self.find_corrections(ypar, dydp, b, tol)
           
            # 4) calculate using the proposed corrections
            egnvls_try = super().getEigenvalues(ypar+x, is_weighted=is_weighted)

            chi2_try, *_ = self.ml.calculate_stats(
                egnvls_try[:, 9], egnvls_try[:,8], egnvls_try[:,11], is_weighted
            )

            self.progress_str.write(
                f'\n\nTOL = {tol}'
                f'\nChi square best = {chisq_best:.8f} from the previous iteration {it-1}'
                f'\nChi square current = {chi2_try:.8f} for the current iteration {it} \n'
            )

            # 5) Try to change the corrections
            is_failed = False

            if chi2_try < chisq_best:
                chisq_best = chi2_try
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

                    egnvls_new = super().getEigenvalues(ypar_try, is_weighted=is_weighted)

                    chi2_new, *_ = self.ml.calculate_stats(
                        egnvls_new[:,9], egnvls_new[:,8], egnvls_new[:,11], is_weighted
                    )

                    if chi2_new < chisq_best:
                        chisq_best = chi2_new
                        ypar = ypar_try

                        self.progress_str.write(
                            f'\nChi square changing corrections = ' \
                            f'{chisq_best:.5f} for iteration {it}'
                        )
                    else:  # chi2new >= chi2best:
                        is_failed = True

                        self.progress_str.write(
                            f'\nChi square final after changing the corrections = ' \
                            f'{chisq_best:.8f} for iteration {it}\n'
                        )

                        break
            # if it isn't get improved then change the tol
            if is_failed:
                tol /= 5.0

            # 6. Save the final parameters
            egnvls = super().getEigenvalues(ypar, is_weighted=is_weighted)

            if self.overwrite:
                Channel.edit_channel_parameters(ypar, self.ml.channels)

                if len(self.ml.couplings) > 0:
                    Coupling.edit_coupling_parameters(ypar, self.ml.couplings)

            _, rms_final, rmsd_final = self.ml.calculate_stats(
                egnvls[:,9], egnvls[:,8], egnvls[:,11], is_weighted
            )

            self.progress_str.write(
                f'\nChi square final/best = {chisq_best:.8f} for iteration {it}\n'
                f'\nRMS final = {rms_final:.8f} cm-1'
                f'\nRMSD final = {rmsd_final:.8f}\n'
            )

            if self.progress_fit:
                print(self.progress_str.getvalue())

            #if rms_fin - rms_init > 0.001:
            #     save_best_parameters(rms_fin)


    def generate_design_matrix_numerical(self, 
        ypar, yfixed, ycal_init, is_weighted):

        """Compute the derivatives matrix by the parameters

        Args:
            ypar (ndarray): fitting parameters
            yfixed (ndarray): fixed/free parameter
            ycal_init (ndarray): the initial calculated values
            is_weighted (bool): apply weighted fit

        Returns:
            ndarray: The matrix with the derivatives
        """

        # TODO: Rewrite to optimize this procedure with numpy
        # TODO: check for some internal python sinc interpolation procedure
        perc = 1.0e-4

        # get changes in the parameters
        dp = np.zeros(ypar.shape[0])
        dydp = np.zeros(shape=(ycal_init.shape[0], ypar.shape[0]), dtype=np.float64)
        
        for prm in range(0, len(ypar)):
            # set all changes to zero
            dp.fill(0.0)

            if yfixed[prm]:
                dp[prm] = abs(ypar[prm]) * perc
                dpar = ypar + dp

                out_data = super().getEigenvalues(dpar, is_weighted=is_weighted)

                ycal = out_data[:,8]
                yvar = out_data[:,11]

                for i in range(0, len(ycal)):
                    dydp[i, prm] = (ycal[i] - ycal_init[i]) / (dp[prm] * yvar[i])
        
        return dydp

    def find_corrections(self, ypar, dydp, b, tol):

        """Find the corrections to the parameters by solving 
        the system of normal equations with SVD technique

        Args:
            ypar (ndarray): the fitting parameters
            dydp (ndarray): The derivatives matrix
            b (ndarray): right hand side vector for solving the normal equations
            tol (float): tolerance value

        Returns:
            ndarray: The suggested corrections to the parameters

        Remarks:
            1. There are two similar procedures for svd fit from numpy 
                and scipy modules In the older versions of the modules there 
                is a difference in the way the default tol value is defined
            2. Tol controls which and how many linear combinations of the 
                parameters will be ignored since the matrix is singular
        """

        #x, _, rank, s = np.linalg.lstsq(dydp, b, tol)
        x, _, rank, s = sp.linalg.lstsq(dydp, b, tol)

        self.progress_str.write(
            f'\nSVD matrix rank = the number significant singular values = {rank}' 
            f'\ntol * max singular value = {tol*s[0]:.3f}\n'
            f'\nThe first {2*rank} significant singular values:\n'
        )
        self.progress_str.write(' '.join(['{0:0.5e}'.format(i) for i in s[0:2*rank]]))
        self.progress_str.write(f'\n\nParameters:\n')
        # TODO: what is some parameter has 0 value?
        astr = np.array2string(np.column_stack((ypar, x, x+ypar, np.absolute(x/ypar)*100)))
        self.progress_str.write(astr)

        return x
   
    def generate_design_matrix_analytical(self, 
        ypar, yfixed, ycal_init, is_weighted):

        if self.ml.sderiv is None:
            raise SystemExit(
                'Spline derivatives should be calculated. Change model to cspline.'
            )

        dydp = np.zeros(shape=(ycal_init.shape[0], ypar.shape[0]), dtype=np.float64)

        evec_selected = np.array([])
        n = 0
        print(self.ml.evecs_index)

        for i in range(0, len(self.ml.evecs_index)):
            evec = np.loadtxt(os.path.join('eigenvectors', f'evector_{i}.dat'))
            n = evec.shape[0]
            evec_selected = np.append(evec_selected, evec[:,self.ml.evecs_index[i]])
        
        evec_selected = evec_selected.reshape(n,-1)

        dydp = np.zeros((evec_selected.shape[1], ypar.shape[0]))

        for prm in range(0, ypar.shape[0]):
            if yfixed[prm]:
                sk = np.diag(self.ml.sderiv[:,prm])
                enr_vec = np.diag(evec_selected.T @ (sk @ evec_selected))
                dydp[:,prm] = enr_vec

        #dydp = dydp * Const.hartree  ## ?????????????????

        np.savetxt('dydp.dta', dydp, fmt='%.6e')
        
        return dydp

class LevMarFit(Fit):

    """Class defining the Levenberg-Marquardt algorithm for nonliner least squares fitting

    Args:
        Fit (object): The inherited class
    """

    def __init__(self, ml):

        super().__init__(ml)

