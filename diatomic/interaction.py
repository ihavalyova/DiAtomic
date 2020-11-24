from math import sqrt as _sqrt
import numpy as np
from utils import C_hartree

__all__ = ['Interaction']


class Interaction:

    def __init__(self, ngrid, nch, ncp, rgrid2):

        self.ngrid = ngrid
        self.nch = nch
        self.rgrid2 = rgrid2
        self.countl = 0
        self.cp_map = {}
        max_rotnj = 100

        # TODO: the upper limit is not correct
        for cp in range(0, ncp+6):
            self.cp_map[cp] = np.zeros((self.nch * self.ngrid, max_rotnj))

    def get_interactions_matrix(self, jrotn, mass, par, channels, couplings,
                                fgrid, dd, Gy2):

        self.jrotn = jrotn
        jjrotn = self.jrotn*(self.jrotn + 1.0)

        pert_matrix = np.zeros((self.nch*self.ngrid)**2).reshape(
            self.nch*self.ngrid, self.nch*self.ngrid
        )

        interact_keys = self.define_interaction_keys()

        for cp in range(0, len(couplings)):
            cprops = zip(
                couplings[cp].interact,
                couplings[cp].coupling,
                couplings[cp].multiplier
            )

            ycs = fgrid[cp*self.ngrid:(cp+1)*self.ngrid]

            for countc, (inter, ctype, m) in enumerate(cprops):
                ch1, ch2 = inter[:2]
                args = self.get_quantum_numbers(channels, ch1, ch2)
                cfunc = interact_keys[ctype](jjrotn, mass, m, par, args)  # Gy2
                row1, row2 = (ch1-1)*self.ngrid, ch1*self.ngrid
                col1, col2 = (ch2-1)*self.ngrid, ch2*self.ngrid

                pert_matrix[row1:row2, col1:col2][dd] += cfunc * ycs
                self.cp_map[cp+countc][row1:row2, self.countl] = cfunc

        self.countl += 1

        return pert_matrix

    def get_quantum_numbers(self, channels, ch1, ch2):

        return {
            'lm1': channels[ch1-1].nlambda,
            'sg1': channels[ch1-1].nsigma,
            'om1': channels[ch1-1].omega,
            'lm2': channels[ch2-1].nlambda,
            'sg2': channels[ch2-1].nsigma,
            'om2': channels[ch2-1].omega,
            's1': (channels[ch1-1].mult-1)/2,
            's2': (channels[ch2-1].mult-1)/2
        }

    def define_interaction_keys(self):

        return {
            'spin-orbit': self.spin_orbit_interaction,
            'lj': self.lj_interaction,
            'lj_parity': self.lj_parity_interaction,
            'sj': self.sj_interaction,
            'sj_parity': self.sj_parity_interaction,
            'sl': self.spin_electornic_interaction,
            'lambdadf': self.lambda_doubling_f_parity,
            'lambdade': self.lambda_doubling_e_parity,
            'lambdad': self.lambda_doubling,
            'ndbobc': self.NDBOBC,
            'dbobc': self.DBOBC,
            'spin-rot': self.spin_rotation_interaction,
            'spin-spin': self.spin_spin_interaction,
        }

    def spin_orbit_interaction(self, jjrotn, mass, m, par, args):

        """
        Calculate diagonal and off-diagonal Spin-Orbit coupling matrix element
        <State1| SO |State1> = multiplier * A(R)
        <State1| SO |State2> = multiplier * alpha(R)

        with selection rules:

        """

        socoef = m * (self.rule_SOdiag(args) or self.rule_SOnondiag(args))

        # return (ycs / C_hartree) * socoef
        return socoef

    def rule_SOdiag(self, args):

        if args['lm1'] == args['lm2'] and \
           args['sg1'] == args['sg2'] and \
           args['om1'] == args['om2']:
            return 1

        return 0

    def rule_SOnondiag(self, args):

        if args['lm1'] != args['lm2'] or \
           args['sg1'] != args['sg2'] or \
           args['om1'] != args['om2']:
            return 1

        return 0

    def lj_interaction(self, jjrotn, mass, m, par, args):

        """
        Calculate the matrix element of L-uncoupling operator
        <State1| LJ |State1> =
        <State1| LJ |State2> =

        with selection rules:
        """

        qexpression = 1.0

        rule_omegaj = self.rule_lj_interaction(args) and \
            _sqrt(jjrotn - (args['om1'] * args['om2']))

        ljcoef = m * qexpression * rule_omegaj

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = -1.0

        return sign * brot * ljcoef

    def rule_lj_interaction(self, args):

        # to check the abs values
        omega_diff = abs(args['om1'] - args['om2'])
        lambda_diff = abs(args['lm1'] - args['lm2'])

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
        sg1 = args['sg1']
        if args['om1'] == 0 and args['sg1'] < 0:
            sg1 = abs(args['sg1'])

        sg2 = args['sg2']
        if args['om2'] == 0 and args['sg2'] < 0:
            sg2 = abs(args['sg2'])

        # if args['sg1'] == args['sg2'] and \
        # args['s1'] == args['s2'] and rule_lj and \

        if sg1 == sg2 and args['s1'] == args['s2'] and rule_lj and \
           (args['om1'] <= self.jrotn and args['om2'] <= self.jrotn+1.0):
            return 1

        return 0

    def lj_parity_interaction(self, jjrotn, mass, m, par, args):

        qexpression = _sqrt(jjrotn + args['om1'] * args['om2'])

        ljcoef = m * qexpression * self.rule_lj_parity(args)

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0  # f-parity
        if par == 1:
            sign = -1.0  # e-parity

        return sign * brot * ljcoef

    def rule_lj_parity(self, args):

        """
        This interaction is allowed only between Sigma-Pi or Sigma-Sigma states
        with even multiplicity with Lambda <= 1 and abs(Omega) = 0.5
        """

        rule_even_mult = (2 * args['sg1'] + 1) % 2 == 0 and \
            (2 * args['sg2'] + 1) % 2 == 0

        if 0 <= args['lm1'] <= 1 and 0 <= args['lm2'] <= 1 and \
           (not args['lm1'] == args['lm2'] == 1) and rule_even_mult and \
           args['om1'] == args['om2'] == 0.5:
            return 1

        return 0

    def sj_interaction(self, jjrotn, mass, m, par, args):

        """
        Calculate the matrix element of S-uncoupling operator
        <State1| SJ |State1> =
        <State1| SJ |State2> =

        with selection rules:
        """

        qexpression = _sqrt(
            args['s1'] * (args['s1'] + 1) - args['sg1'] * args['sg2']
        )

        # TODO: check this!
        rule_omegaj = \
            self.rule_sj_interaction(args) and \
            _sqrt(jjrotn - (args['om1'] * args['om2']))

        sjcoef = m * qexpression * rule_omegaj

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = -1.0

        return sign * brot * sjcoef

    def rule_sj_interaction(self, args):

        omega_diff = args['om1'] - args['om2']
        sigma_diff = args['sg1'] - args['sg2']
        rule_sj = \
            (omega_diff == sigma_diff == 1.0) or \
            (omega_diff == sigma_diff == -1.0)

        if args['lm1'] == args['lm2'] and args['s1'] == args['s2'] and \
           (args['om1'] <= self.jrotn and args['om2'] <= self.jrotn+1.0) and \
           rule_sj:
            return 1

        return 0

    def sj_parity_interaction(self, jjrotn, mass, m, par, args):
        """
            only between Sigma states
        """

        qexpression = _sqrt(jjrotn + args['om1'] * args['om2'])

        sjcoef = m * qexpression * self.rule_sj_parity(args)

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0  # f-parity
        if par == 1:
            sign = -1.0  # e-parity

        return sign * brot * sjcoef

    def rule_sj_parity(self, args):

        if args['lm1'] == args['lm2'] == 0 and args['s1'] == args['s2']:
            return 1

        return 0

    def spin_electornic_interaction(self, jjrotn, mass, m, par, args):
        """
        Calculate the matrix element of spin-electronic operator
        <State1| LS |State1> =
        <State1| LS |State2> =

        with selection rules:
        """

        qexpression = _sqrt(
            (args['s1'] * (args['s1'] + 1)) - args['sg1'] * args['sg2']
        )

        slcoef = m * qexpression * self.rule_spin_electronic(args)

        brot = 1.0 / (2.0 * mass * self.rgrid2)

        sign = 1.0

        return sign * brot * slcoef

    def rule_spin_electronic(self, args):

        lambda_diff = args['lm1'] - args['lm2']
        sigma_diff = args['sg2'] - args['sg1']

        rule_ls = \
            (lambda_diff == sigma_diff == 1.0) or \
            (lambda_diff == sigma_diff == -1.0)

        if args['om1'] == args['om2'] and args['s1'] == args['s2'] and rule_ls:
            return 1

        return 0

    def lambda_doubling(self, jjrotn, mass, m, par, args):

        ldcoef = jjrotn / (2.0 * mass * self.rgrid2)**2

        return m * ldcoef

    def lambda_doubling_e_parity(self, jjrotn, mass, m, par, args):

        if par == 1:
            return self.lambda_doubling(jjrotn, mass, m, par, args)

        return np.zeros(self.rgrid2.shape[0])

    def lambda_doubling_f_parity(self, jjrotn, mass, m, par, args):

        if par == 0:
            return self.lambda_doubling(jjrotn, mass, m, par, args)

        return np.zeros(self.rgrid2.shape[0])

    def NDBOBC(self, jjrotn, mass, m, par, args):

        bocoef = jjrotn / (2.0 * mass * self.rgrid2)

        return m * bocoef * C_hartree

    def DBOBC(self, jjrotn, mass, m, par, args):

        return (1.0 / C_hartree) * m

    def spin_rotation_interaction(self, jjrotn, mass, m, par, args):

        if self.rule_spin_rot_diag(args):

            qexpression = args['sg1']**2 - args['s1']*(args['s2'] + 1.0)

        elif self.rule_spin_rot_nondiag(args):
            ss1 = args['s1'] * (args['s1'] + 1)
            qexpression = \
                _sqrt(jjrotn - (args['om1'] * args['om2'])) * \
                _sqrt(ss1 - args['sg1'] * args['sg2'])

        srcoef = m * qexpression

        sign = 1.0

        return sign * srcoef

    def rule_spin_rot_diag(self, args):

        if args['om1'] == args['om2'] and \
           args['s1'] == args['s2'] and \
           args['sg1'] == args['sg2'] and \
           args['lm1'] == args['lm2']:
            return 1

        return 0

    def rule_spin_rot_nondiag(self, args):

        self.rule_sj_interaction(args)

    def spin_spin_interaction(self, jjrotn, mass, m, par, args):

        if self.rule_spin_spin_diag(args):
            qexpression = 3 * args['sg1']**2 - args['s1']*(args['s2'] + 1.0)

        elif self.rule_spin_spin_nondiag(args):
            qexpression = 1.0

        sscoef = m * qexpression

        sign = 1.0

        return sign * sscoef

    def rule_spin_spin_nondiag(self, args):

        lambda_diff = args['lm1'] - args['lm2']
        sigma_diff = args['sg2'] - args['sg1']

        # Delta Sigma = +/- 1
        rule_ss1_1so = \
            (lambda_diff == -1.0 * sigma_diff == 1.0) or \
            (lambda_diff == -1.0 * sigma_diff == -1.0)

        # Delta Sigma = +/- 2
        rule_ss2_2so = \
            (lambda_diff == -1.0 * sigma_diff == 2.0) or \
            (lambda_diff == -1.0 * sigma_diff == -2.0)

        rule_ss1 = \
            (args['s1'] - args['s2']) == 1.0 or \
            (args['s1'] - args['s2']) == -1.0 and \
            rule_ss1_1so

        rule_ss2 = \
            (args['s1'] - args['s2']) == 2.0 or \
            (args['s1'] - args['s2']) == -2.0 and \
            rule_ss2_2so

        if rule_ss1 or rule_ss2:
            return 1

        return 0

    def rule_spin_spin_diag(self, args):

        if args['om1'] == args['om2'] and \
           args['s1'] == args['s2'] and \
           args['sg1'] == args['sg2']:
            return 1

        return 0
