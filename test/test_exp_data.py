from diatom import MoleculeData
import unittest
import io


class TestExpData(unittest.TestCase):

    def test_existing_exp_file(self):

        mdata = MoleculeData()
        self.assertRaises(
            SystemExit, mdata.set_exp_data, 'some_file', markers=[1]
        )

    def test_exp_data_shape_large(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5
        exp_data_string = u"""8 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        4  0  5.5  107.6598777  1  2.993E-03  1  5
        5  0  6.5  159.4241657  1  3.105E-03  1  5
        6  0  7.5  219.0977590  1  3.052E-03  1  5
        7  0  8.5  286.6573545  1  3.059E-03  1  5
        8  0  9.5  362.0793175  1  3.108E-03  1  5
        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream, markers=[1])

        self.assertEqual(mdata.exp_data.shape[0], 8)

    def test_exp_data_shape2(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5
        exp_data_string = u"""2 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream, markers=[1])

        self.assertEqual(mdata.exp_data.shape[0], 2)

    def test_exp_data_shape1(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5,
        exp_data_string = u"""1 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream, markers=[1])

        self.assertEqual(mdata.exp_data.shape[0], 1)

    def test_exp_data_markers_None(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        exp_data_string = u"""3 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream)

        self.assertEqual(mdata.exp_data.shape[0], 3)

    def test_exp_data_markers_wrong_value(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        exp_data_string = u"""3 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        """

        data_stream = io.StringIO(exp_data_string)

        self.assertRaises(
            SystemExit, mdata.set_exp_data, data_stream, markers=[2]
        )

    def test_exp_data_jnumbers_wrong_values(self):

        mdata = MoleculeData()
        mdata.jvalues = 8.5, 9.5, 6.5
        exp_data_string = u"""3 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        """

        data_stream = io.StringIO(exp_data_string)

        self.assertRaises(
            SystemExit, mdata.set_exp_data, data_stream, markers=[1]
        )

    def test_exp_data_both_parities(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        mdata.parities = 0, 1
        exp_data_string = u"""6 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        4  0  2.5  0.000000000  0  5.0e-03    2  5
        5  0  3.5  27.92907799  0  2.554E-03  2  5
        6  0  4.5  63.82044637  0  2.589E-03  2  5
        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream, markers=[1, 2])

        self.assertTrue(mdata.exp_data.shape[0] == 6)

    def test_exp_data_missing_first_line(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        exp_data_string = u"""
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        """

        data_stream = io.StringIO(exp_data_string)

        self.assertRaises(
            SystemExit, mdata.set_exp_data, data_stream, markers=[1]
        )

    def test_exp_data_missing_column(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        exp_data_string = u"""3 \n
        1  0  2.5  0.000000000  1  5.0e-03    1
        2  0  3.5  27.92907799  1  2.554E-03  1
        3  0  4.5  63.82044637  1  2.589E-03  1
        """

        data_stream = io.StringIO(exp_data_string)

        self.assertRaises(
            SystemExit, mdata.set_exp_data, data_stream, markers=[1]
        )

    def test_exp_data_wrong_number_of_lines(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        exp_data_string = u"""15 \n
        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream)

        self.assertEqual(mdata.exp_data.shape[0], 3)

    def test_exp_data_empty_lines(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        exp_data_string = u"""

        3 \n

        1  0  2.5  0.000000000  1  5.0e-03    1  5
        2  0  3.5  27.92907799  1  2.554E-03  1  5

        3  0  4.5  63.82044637  1  2.589E-03  1  5

        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream)

        self.assertEqual(mdata.exp_data.shape[0], 3)

    def test_exp_data_comments(self):

        mdata = MoleculeData()
        mdata.jvalues = 2.5, 3.5, 4.5
        exp_data_string = u"""3 \n                  # comment
        # comment
        # comment
        1  0  2.5  0.000000000  1  5.0e-03    1  5  # comment
        # comment
        2  0  3.5  27.92907799  1  2.554E-03  1  5
        3  0  4.5  63.82044637  1  2.589E-03  1  5
        #4  0  4.5  64.82044637  1  2.589E-03  1  5   # skip commented lines
        #5  0  4.5  65.82044637  1  2.589E-03  1  5   # skip commented lines
        # comment
        # comment
        """

        data_stream = io.StringIO(exp_data_string)
        mdata.set_exp_data(data_stream)

        self.assertEqual(mdata.exp_data.shape[0], 3)
