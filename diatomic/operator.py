import numpy as np

from math import sqrt as _sqrt
from scipy.interpolate import CubicSpline as _CubicSpline

from .interpolator import CSpline
from .utils import Utils as _utils

__all__ = ['Operator', 'SO', 'LJ', 'SJ', 'LS', 'LD', 'Dipole']


class Operator(object):

    # TODO: introduce a units parameter for each operator to
    # convert all parameters in ypar to the same units.

    def __init__(self, objs=None, pair_states=None, label='', model='pointwise',
                 rotc=0.0, multiplier=1, custom_func=None, shift_by=0.0):

        diatomic, grid, states = objs[0], objs[1], objs[2]

        # get grid data
        self.ngrid = grid.ngrid
        self.rmin = grid.rmin
        self.rmax = grid.rmax
        self.rgrid = grid.rgrid
        self.rgrid2 = np.square(self.rgrid)
        self.rstep = grid.rstep
        self.solver = grid.solver.lower()

        self.molecule = diatomic.molecule
        self.mol_name = ''.join(
            filter(lambda x: not x.isdigit(), self.molecule[0]))

        # get diatomic data
        self.masses = diatomic.reduced_masses or diatomic.masses
        self.niso = diatomic.niso
        self.refj = diatomic.referencej
        self.ref_enr = diatomic.ref_enr
        self.exp_data = diatomic.exp_data
        self.exp_file = diatomic.exp_file
        self.wavens_data = diatomic.wavens_data
        self.params_by_labels = diatomic.params_by_labels
        self.labels_inds = diatomic.labels_inds
        self.params = diatomic.params
        self.fixed = diatomic.fixed
        self.fname_data_params = diatomic.fname_data_params

        # set the local parameters
        self.states = states
        self.nch = len(self.states)
        self.msize = self.nch * self.ngrid
        self.pair_states = pair_states
        self.label = label
        self.model = model.lower()
        self.rot_correction = rotc
        self.multiplier = multiplier
        self.custom_func = custom_func
        self.shift_by = shift_by
        self.xunits = 1.0 / _utils.C_bohr
        self.yunits = 1.0 / _utils.C_hartree
        self.ygrid = np.zeros(self.ngrid)
        self.matrix = np.zeros((self.ngrid, self.ngrid))
        self.dd = np.diag_indices(self.ngrid)

    def set_states(self):
        # the indices of the states
        self.istate1 = self.pair_states[0]
        self.istate2 = self.pair_states[1]

        # the state objects
        for state in self.states:
            if state.nstate == self.istate1:
                self.state1 = state
            if state.nstate == self.istate2:
                self.state2 = state


class SO(Operator):

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)

        self.set_states()

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        self.matrix[self.dd] = self.spin_orbit_interaction() * self.ygrid

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_coupling_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_coupling_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_coupling_on_grid(ypar)

    def _calculate_pointwise_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # The parameters in ypar are in cm-1 units
        ypnts = ypar[self.sind:self.eind] * self.yunits  # convert to au units
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk = cs(self.rgrid, return_deriv=True)

    def _calculate_custom_coupling_on_grid(self, ypar):

        xpnts = self.rgrid
        # The parameters in ypar are in au units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        self.ygrid = self.custom_func(xpnts, ypnts)

    def spin_orbit_interaction(self):
        """Calculate diagonal and off-diagonal spin-orbit matrix element as:
            ``<State1| SO |State1>`` = multiplier * A(R)
            ``<State1| SO |State2>`` = multiplier * alpha(R)

        Returns:
            array: the computed matrix elements
        """

        soc = self.multiplier * (self.rule_so_diag() or self.rule_so_nondiag())

        # return (ycs / C_hartree) * socoef
        return soc

    def rule_so_diag(self):

        if self.state1._lambda == self.state2._lambda and \
           self.state1.sigma == self.state2.sigma and \
           self.state1.omega == self.state2.omega:
            return 1

        return 0

    def rule_so_nondiag(self):

        if self.state1._lambda != self.state2._lambda or \
           self.state1.sigma != self.state2.sigma or \
           self.state1.omega != self.state2.omega:
            return 1

        return 0


class LJ(Operator):

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)

        self.set_states()

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        mass = self.masses[iso-1]
        jjrotn = jrotn*(jrotn + 1.0)

        # this check if not complete
        if self.state1.omega == self.state2.omega == 0.5:
            self.matrix[self.dd] = \
                self._lj_parity_interaction(jjrotn, mass, par) * self.ygrid
        else:
            self.matrix[self.dd] = \
                self._lj_interaction(jrotn, jjrotn, mass) * self.ygrid

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_coupling_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_coupling_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_coupling_on_grid(ypar)

    def _calculate_pointwise_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # The parameters in ypar are in cm-1 units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk = cs(self.rgrid, return_deriv=True)

    def _calculate_custom_coupling_on_grid(self, ypar):

        xpnts = self.rgrid
        # The parameters in ypar are in au units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        self.ygrid = self.custom_func(xpnts, ypnts)

    def _lj_interaction(self, jrotn, jjrotn, mass):
        """Calculate the matrix elements of LJ operator as:
           ``<State1| LJ |State1>`` =
           ``<State1| LJ |State2>`` =

        Args:
            jrotn (float): the value of J
            mass (float): the value of mass

        Returns:
            array: the computed matrix elements
        """

        qexpression = 1.0
        rule_omegaj = self._rule_lj_interaction(jrotn) and \
            _sqrt(jjrotn - (self.state1.omega * self.state2.omega))

        lj = self.multiplier * qexpression * rule_omegaj
        brot = 1.0 / (2.0 * mass * self.rgrid2)
        sign = -1.0

        return sign * brot * lj

    def _rule_lj_interaction(self, jrotn):

        # to check the abs values
        omega_diff = abs(self.state1.omega - self.state2.omega)
        lambda_diff = abs(self.state1._lambda - self.state2._lambda)

        rule_lj = \
            ((omega_diff == lambda_diff) == 1.0) or \
            ((omega_diff == lambda_diff) == -1.0)

        # rule_lj = \
        #   ((args['om1'] - args['om2']) == \
        #   (args['lm1'] - args['lm2']) == 1.0) or \
        #    ((args['om1'] - args['om2']) == \
        #   (args['lm1'] - args['lm2']) == -1.0)
        # print('LJ', rule_lj)

        # # to check!!!
        sg1 = self.state1.sigma
        if self.state1.omega == 0 and self.state1.sigma < 0:
            sg1 = abs(self.state1.sigma)

        sg2 = self.state2.sigma
        if self.state2.omega == 0 and self.state2.sigma < 0:
            sg2 = abs(self.state2.sigma)

        # if args['sg1'] == args['sg2'] and \
        # args['s1'] == args['s2'] and rule_lj and \

        if sg1 == sg2 and self.state1.spin == self.state2.spin and rule_lj and \
           (self.state1.omega <= jrotn and self.state2.omega <= jrotn+1.0):
            return 1

        return 0

    def _lj_parity_interaction(self, jjrotn, mass, par):

        qexpression = _sqrt(jjrotn + self.state1.omega * self.state1.omega)
        ljcoef = self.multiplier * qexpression * self._rule_lj_parity()
        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0  # f-parity
        if par == 1:
            sign = -1.0  # e-parity

        return sign * brot * ljcoef

    def _rule_lj_parity(self):
        """This kind of interaction is allowed only between Sigma-Pi or
        Sigma-Sigma states with even multiplicity with Lambda <= 1 and
        abs(Omega) = 0.5
        """

        rule_even_mult = (2 * self.state1.sigma + 1) % 2 == 0 and \
            (2 * self.state2.sigma + 1) % 2 == 0

        if 0 <= self.state1._lambda <= 1 and \
           0 <= self.state2._lambda <= 1 and \
           (not self.state1._lambda == self.state2._lambda == 1) and \
           rule_even_mult and self.state1.omega == self.state2.omega == 0.5:
            return 1

        return 0


class SJ(Operator):

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)

        self.set_states()

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        mass = self.masses[iso-1]
        jjrotn = jrotn*(jrotn + 1.0)

        if self.istate1 != self.istate2:
            self.matrix[self.dd] = \
                self._sj_interaction(jrotn, jjrotn, mass) * self.ygrid
        else:
            self.matrix[self.dd] = \
                self._sj_parity_interaction(jjrotn, mass, par) * self.ygrid

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_coupling_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_coupling_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_coupling_on_grid(ypar)

    def _calculate_pointwise_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # The parameters in ypar are in cm-1 units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk = cs(self.rgrid, return_deriv=True)

    def _calculate_custom_coupling_on_grid(self, ypar):

        xpnts = self.rgrid
        # The parameters in ypar are in au units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        self.ygrid = self.custom_func(xpnts, ypnts)

    def _sj_interaction(self, jrotn, jjrotn, mass):
        """Calculate the matrix elements of SJ operator which are:
            ``<State1| SJ |State1>`` =
            ``<State1| SJ |State2>`` =

        Args:
            jrotn (float): the value of J
            mass (float): the value of mass

        Returns:
            array: the computed matrix elements
        """
        ss1 = self.state1.spin * (self.state1.spin + 1)
        qexpression = _sqrt(ss1 - self.state1.sigma * self.state2.sigma)

        # TODO: check this!
        rule_omegaj = \
            self._rule_sj_interaction(jrotn) and \
            _sqrt(jjrotn - (self.state1.omega * self.state2.omega))

        sjcoef = self.multiplier * qexpression * rule_omegaj
        brot = 1.0 / (2.0 * mass * self.rgrid2)
        sign = -1.0

        return sign * brot * sjcoef

    def _rule_sj_interaction(self, jrotn):

        omega_diff = self.state1.omega - self.state2.omega
        sigma_diff = self.state1.sigma - self.state2.sigma
        rule_sj = \
            (omega_diff == sigma_diff == 1.0) or \
            (omega_diff == sigma_diff == -1.0)

        if self.state1._lambda == self.state2._lambda and \
           self.state1.spin == self.state2.spin and \
           (self.state1.omega <= jrotn and self.state2.omega <= jrotn+1.0) and \
           rule_sj:
            return 1

        return 0

    def _sj_parity_interaction(self, jjrotn, mass, par):
        """only between Sigma states"""

        qexpression = _sqrt(jjrotn + self.state1.omega * self.state2.omega)
        sjcoef = self.multiplier * qexpression * self._rule_sj_parity()
        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0  # f-parity
        if par == 1:
            sign = -1.0  # e-parity

        return sign * brot * sjcoef

    def _rule_sj_parity(self):

        if self.state1._lambda == self.state2._lambda == 0 and \
           self.state1.spin == self.state2.spin:
            return 1

        return 0


class LS(Operator):

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)

        self.set_states()

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        mass = self.masses[iso-1]

        self.matrix[self.dd] = \
            self._spin_electornic_interaction(mass) * self.ygrid

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_coupling_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_coupling_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_coupling_on_grid(ypar)

    def _calculate_pointwise_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # The parameters in ypar are in cm-1 units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk = cs(self.rgrid, return_deriv=True)

    def _calculate_custom_coupling_on_grid(self, ypar):

        xpnts = self.rgrid
        # The parameters in ypar are in au units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        self.ygrid = self.custom_func(xpnts, ypnts)

    def _spin_electornic_interaction(self, mass):
        """Calculate the matrix elements of SL operator which are:
            ``<State1| SL |State1>`` =
            ``<State1| SL |State2>`` =

        Args:
            mass (float): the value of the mass

        Returns:
            array: the computed matrix elements
        """

        ss1 = self.state1.spin * (self.state2.spin + 1)
        qexpression = _sqrt(ss1 - self.state1.sigma * self.state2.sigma)
        sl = self.multiplier * qexpression * self._rule_spin_electronic()
        brot = 1.0 / (2.0 * mass * self.rgrid2)
        sign = 1.0

        return sign * brot * sl

    def _rule_spin_electronic(self):

        lambda_diff = self.state1._lambda - self.state2._lambda
        sigma_diff = self.state2.sigma - self.state1.sigma

        rule_ls = \
            (lambda_diff == sigma_diff == 1.0) or \
            (lambda_diff == sigma_diff == -1.0)

        if self.state1.omega == self.state2.omega and \
           self.state1.spin == self.state2.spin and \
           rule_ls:
            return 1

        return 0


class LD(Operator):

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)
        self.set_states()

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        mass = self.masses[iso-1]
        jjrotn = jrotn * (jrotn + 1.0)

        self.matrix[self.dd] = \
            self._lambda_doubling(jjrotn, mass, par) * self.ygrid

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_coupling_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_coupling_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_coupling_on_grid(ypar)

    def _calculate_pointwise_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # The parameters in ypar are in cm-1 units
        ypnts = ypar[self.sind:self.eind] / self.yunits
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk = cs(self.rgrid, return_deriv=True)

    def _calculate_custom_coupling_on_grid(self, ypar):

        xpnts = self.rgrid
        # The parameters in ypar are in au units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        self.ygrid = self.custom_func(xpnts, ypnts)

    def _lambda_doubling(self, jjrotn, mass, par):

        ldcoef = jjrotn / (2.0 * mass * self.rgrid2)**2

        return self.multiplier * ldcoef

    def _lambda_doubling_e_parity(self, jjrotn, mass, par):

        if par == 1:
            return self._lambda_doubling(jjrotn, mass, par)

        return np.zeros(self.rgrid2.shape[0])

    def _lambda_doubling_f_parity(self, jjrotn, mass, par):

        if par == 0:
            return self._lambda_doubling(jjrotn, mass, par)

        return np.zeros(self.rgrid2.shape[0])


class SR(Operator):

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)
        self.set_states()

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        # mass = self.masses[iso-1]
        jjrotn = jrotn * (jrotn + 1.0)

        self.matrix[self.dd] = \
            self._spin_rotation_interaction(jrotn, jjrotn) * self.ygrid

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_coupling_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_coupling_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_coupling_on_grid(ypar)

    def _calculate_pointwise_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # The parameters in ypar are in cm-1 units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk = cs(self.rgrid, return_deriv=True)

    def _calculate_custom_coupling_on_grid(self, ypar):

        xpnts = self.rgrid
        # The parameters in ypar are in au units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        self.ygrid = self.custom_func(xpnts, ypnts)

    def _spin_rotation_interaction(self, jrotn, jjrotn):

        if self._rule_spin_rot_diag():
            qexpr = self.state1.sigma**2-self.state1.spin*(self.state2.spin+1)
        elif self._rule_spin_rot_nondiag(jrotn):
            ss1 = self.state1.spin * (self.state1.spin + 1)
            qexpr = \
                _sqrt(jjrotn - (self.state1.omega * self.state2.omega)) * \
                _sqrt(ss1 - self.state1.sigma * self.state2.sigma)

        srcoef = self.multiplier * qexpr
        sign = 1.0

        return sign * srcoef

    def _rule_spin_rot_diag(self):

        if self.state1.omega == self.state2.omega and \
           self.state1.spin == self.state2.spin and \
           self.state1.sigma == self.state2.sigma and \
           self.state1._lambda == self.state2._lambda:
            return 1

        return 0

    def _rule_spin_rot_nondiag(self, jrotn):
        self._rule_sj_interaction(jrotn)

    def _rule_sj_interaction(self, jrotn):

        omega_diff = self.state1.omega - self.state2.omega
        sigma_diff = self.state1.sigma - self.state2.sigma
        rule_sj = \
            (omega_diff == sigma_diff == 1.0) or \
            (omega_diff == sigma_diff == -1.0)

        if self.state1._lambda == self.state2._lambda and \
           self.state1.spin == self.state2.spin and \
           (self.state1.omega <= jrotn and self.state2.omega <= jrotn+1.0) and \
           rule_sj:
            return 1

        return 0


class SS(Operator):

    def __init__(self, objs, pair_states, label, model='pointwise', rotc=0.0,
                 multiplier=1, custom_func=None, shift_by=0.0):

        Operator.__init__(self, objs, pair_states, label, model=model,
                          rotc=rotc, multiplier=multiplier,
                          custom_func=custom_func, shift_by=shift_by)
        self.set_states()

        # set the radial parameters associated with this operator
        self.init_params = self.params_by_labels[self.label]
        self.sind, self.eind = self.labels_inds[self.label]

        nrows = self.init_params.shape[0]
        self.xunits = np.full(nrows, self.xunits)
        self.yunits = np.full(nrows, self.yunits)
        self.xpoints = self.init_params[:, 0]
        self.ypoints = self.init_params[:, 1]
        self.fixed = self.init_params[:, 2]

        # self.calculate_radial_function_on_grid(self.params)

    def calculate_matrix_elements(self, jrotn, par, iso):

        self.matrix[self.dd] = self._spin_spin_interaction() * self.ygrid

        return self.matrix

    def calculate_radial_function_on_grid(self, ypar):

        if self.model in ['pw', 'pointwise']:
            self._calculate_pointwise_coupling_on_grid(ypar)

        if self.model == 'cspline':
            self._calculate_cspline_coupling_on_grid(ypar)

        if self.model == 'custom':
            self._calculate_custom_coupling_on_grid(ypar)

    def _calculate_pointwise_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        # The parameters in ypar are in cm-1 units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        cs = _CubicSpline(xpnts, ypnts, bc_type='natural', extrapolate=True)
        self.ygrid = cs(self.rgrid)

    def _calculate_cspline_coupling_on_grid(self, ypar):

        xpnts = self.xpoints * self.xunits
        ypnts = ypar[self.sind:self.eind] * self.yunits
        cs = CSpline(xpnts, ypnts)
        self.ygrid, sk = cs(self.rgrid, return_deriv=True)

    def _calculate_custom_coupling_on_grid(self, ypar):

        xpnts = self.rgrid
        # The parameters in ypar are in au units
        ypnts = ypar[self.sind:self.eind]  # * self.yunits
        self.ygrid = self.custom_func(xpnts, ypnts)

    def _spin_spin_interaction(self):

        if self._rule_spin_spin_diag():
            ss1 = self.state1.spin * (self.state2.spin + 1.0)
            qexpression = 3 * self.state1.sigma**2 - ss1

        elif self._rule_spin_spin_nondiag():
            qexpression = 1.0

        sscoef = self.multiplier * qexpression
        sign = 1.0

        return sign * sscoef

    def _rule_spin_spin_nondiag(self):

        lambda_diff = self.state1._lambda - self.state2._lambda
        sigma_diff = self.state2.sigma - self.state1.sigma

        # Delta Sigma = +/- 1
        rule_ss1_1so = \
            (lambda_diff == -1.0 * sigma_diff == 1.0) or \
            (lambda_diff == -1.0 * sigma_diff == -1.0)

        # Delta Sigma = +/- 2
        rule_ss2_2so = \
            (lambda_diff == -1.0 * sigma_diff == 2.0) or \
            (lambda_diff == -1.0 * sigma_diff == -2.0)

        rule_ss1 = \
            (self.state1.spin - self.state2.spin) == 1.0 or \
            (self.state1.spin - self.state2.spin) == -1.0 and \
            rule_ss1_1so

        rule_ss2 = \
            (self.state1.spin - self.state2.spin) == 2.0 or \
            (self.state1.spin - self.state2.spin) == -2.0 and \
            rule_ss2_2so

        if rule_ss1 or rule_ss2:
            return 1

        return 0

    def _rule_spin_spin_diag(self):

        if self.state1.omega == self.state2.omega and \
           self.state1.spin == self.state2.spin and \
           self.state1.sigma == self.state2.sigma:
            return 1

        return 0


class Dipole(Operator):

    def __init__(self, first_channel, second_channel, dm_file):

        self.ch1 = first_channel
        self.ch2 = second_channel

    def write_pointwise_data(self, file_name, x, y, z):

        data = np.column_stack((x, y, z))
        fmt = ['%20.10f', '%25.10f', '%6d']
        np.savetxt(file_name, data, fmt=fmt)
