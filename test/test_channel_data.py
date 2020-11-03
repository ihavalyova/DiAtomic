from diatom import Channel
from diatom import Const
import unittest
import io
import warnings
import numpy as np


class TestChannelData(unittest.TestCase):

    def create_channel(self, filep, model):

        return Channel(
            filep=filep,
            model=model,
            nlambda=0,
            sigma=0.5,
            multiplicity=2,
            rot_correction=0.0
        )

    def test_channel_constructor(self):

        filep = '0sigma_s.pot'
        model = 'pointwise'
        ch = self.create_channel(filep, model)

        self.assertEqual(ch.filep, filep)
        self.assertEqual(ch.model, model)
        self.assertEqual(ch.nlambda, 0)
        self.assertEqual(ch.nsigma, 0.5)
        self.assertEqual(ch.mult, 2)
        self.assertEqual(ch.rot_correction, .0)

    def test_channel_invalid_file(self):

        filep = 'invalid.pot'
        model = 'pointwise'
        ch = self.create_channel(filep, model)

        self.assertRaises(
            OSError, Channel.set_channel_parameters, [ch]
        )

    def test_channel_invalid_model(self):

        filep_str = """6 \n
        1.116847826   13498.734441836  0
        1.214673913    7882.831966143  0
        1.508152173    1627.900567243  0
        1.630434782    2353.243138193  0
        1.752717391    3933.653801380  1
        2.085937500   10009.411380241  1
        """
        model = 'invalid'
        filep = io.StringIO(filep_str)
        ch = self.create_channel(filep, model)

        self.assertRaises(
            SystemExit, Channel.set_channel_parameters, [ch]
        )

    def test_channel_pointwise_empty_file(self):

        filep_str = "\n"
        model = 'pointwise'
        filep = io.StringIO(filep_str)
        ch = self.create_channel(filep, model)
        warnings.simplefilter("ignore")

        self.assertRaises(
            IndexError, Channel.set_channel_parameters, [ch]
        )

    def test_channel_pointwise_syntax_file1(self):

        filep_str = """8 \n
        1.116847826   13498.734441836  0
        1.214673913    7882.831966143  0
        1.385869565    2334.929989191  0
        1.508152173    1627.900567243  0
        1.630434782    2353.243138193  0
        1.752717391    3933.653801380  1
        1.972826087    7712.887614177  1
        2.085937500   10009.411380241  1
        """

        model = 'pointwise'
        filep = io.StringIO(filep_str)
        ch = self.create_channel(filep, model)
        Channel.set_channel_parameters([ch])

        self.assertEqual(ch.R.shape[0], 8)
        self.assertEqual(ch.U.shape[0], 8)
        self.assertEqual(ch.fixedU.shape[0], 8)
        self.assertEqual(ch.npnts, 8)

    def test_channel_pointwise_syntax_file2(self):

        filep_str = """8 \n
        1.116847826   13498.734441836
        1.214673913    7882.831966143
        1.385869565    2334.929989191
        1.508152173    1627.900567243
        1.630434782    2353.243138193
        1.752717391    3933.653801380
        1.972826087    7712.887614177
        2.085937500   10009.411380241
        """

        model = 'pointwise'
        filep = io.StringIO(filep_str)
        ch = self.create_channel(filep, model)
        Channel.set_channel_parameters([ch])

        self.assertEqual(ch.R.shape[0], 8)
        self.assertEqual(ch.U.shape[0], 8)
        self.assertEqual(ch.fixedU.shape[0], 8)
        self.assertEqual(ch.npnts, 8)

    def test_channel_pointwise_syntax_file3(self):

        filep_str = """8 \n
        1.116847826
        1.214673913    7882.831966143
        1.385869565    2334.929989191
        1.508152173    1627.900567243
        1.630434782    2353.243138193
        1.752717391    3933.653801380
        1.972826087    7712.887614177
        2.085937500   10009.411380241
        """

        model = 'pointwise'
        filep = io.StringIO(filep_str)
        ch = self.create_channel(filep, model)

        self.assertRaises(
            ValueError, Channel.set_channel_parameters, [ch]
        )

    def test_channel_pointwise_file_comments(self):

        filep_str = """8 # comment \n  # comment
        # comment
        1.116847826   13498.734441836  # comment
        1.214673913    7882.831966143
        # comment
        1.385869565    2334.929989191
        1.508152173    1627.900567243
        1.630434782    2353.243138193  # comment
        # comment
        1.752717391    3933.653801380
        1.972826087    7712.887614177
        #2.085937500   10009.411380241 # comment

        # comment
        """

        model = 'pointwise'
        filep = io.StringIO(filep_str)
        ch = self.create_channel(filep, model)
        Channel.set_channel_parameters([ch])

        self.assertEqual(ch.R.shape[0], 7)
        self.assertEqual(ch.U.shape[0], 7)
        self.assertEqual(ch.fixedU.shape[0], 7)
        self.assertEqual(ch.npnts, 7)

    def test_channel_pointwise_parameters_value(self):

        filep = '0sigma_s.pot'
        model = 'pointwise'
        ch = self.create_channel(filep, model)
        data = np.loadtxt(filep, skiprows=1)

        Channel.set_channel_parameters([ch])

        self.assertTrue(np.allclose(ch.R * Const.bohr, data[:, 0]))
        self.assertTrue(np.allclose(ch.U * Const.hartree, data[:, 1]))
        self.assertTrue(np.array_equal(ch.fixedU, data[:, 2]))
