from diatom import MoleculeData
import unittest
import numpy as np


class TestMoleculeData(unittest.TestCase):

    # def test_reduced_mass_type(self):
    #     mdata = MoleculeData()
    #     mdata._calculate_molecule_reduced_mass(1)

    def test_molecule_symbol_fail(self):

        mdata = MoleculeData()

        # pass
        self.assert_is_instance_of_float(mdata, '1H 1H')
        self.assert_is_instance_of_float(mdata, '1 H 1 H')
        self.assert_is_instance_of_float(mdata, '1 H 1H')
        self.assert_is_instance_of_float(mdata, '1H 1 H')
        self.assert_is_instance_of_float(mdata, ' 1H1H')
        self.assert_is_instance_of_float(mdata, '1H1H ')

        # fail
        self.assert_raised_error(mdata, 'H1H1')
        self.assert_raised_error(mdata, '1HH')
        self.assert_raised_error(mdata, 'H1H1H')

        self.assertTrue(
            mdata._calculate_molecule_reduced_mass('1H1H') < 1.0
        )

        self.assertRaises(
            AttributeError, mdata._calculate_molecule_reduced_mass, 1
        )

    def assert_raised_error(self, mdata, symbol):

        self.assertRaises(
            SystemExit, mdata._calculate_molecule_reduced_mass, symbol
        )

    def assert_is_instance_of_float(self, mdata, symbol):

        self.assertIsInstance(
            mdata._calculate_molecule_reduced_mass(symbol), float
        )

    def test_molecule_property(self):

        mdata = MoleculeData()

        # string should not work
        try:
            mdata.molecule = '58Ni1H'
        except SystemExit:
            pass

        # tuple or list of strings should work
        mdata.molecule = '58Ni1H',
        self.assertTrue(len(mdata.reduced_masses) == 1)

        mdata.molecule = ['58Ni1H']
        self.assertTrue(len(mdata.reduced_masses) == 1)

    def test_masses_property(self):

        mdata = MoleculeData()

        # number or string should not work
        try:
            mdata.masses = 1.0
        except TypeError:
            pass

        # tuple or list of numbers should work
        mdata.masses = 1.0,
        self.assertEqual(len(mdata.imasses), 1)

        mdata.masses = [1.0]
        self.assertEqual(len(mdata.imasses), 1)

    def test_jrange_and_jvalues_property(self):

        mdata = MoleculeData()

        mdata.jrange = 1, 5
        self.assertCountEqual(mdata.jqnumbers, np.array([1, 2, 3, 4, 5]))
        mdata.jqnumbers = np.array([], dtype=np.float64)

        mdata.jrange = 1, 1
        self.assertCountEqual(mdata.jqnumbers, np.array([1]))
        mdata.jqnumbers = np.array([], dtype=np.float64)

        mdata.jrange = 1, 0
        self.assertEqual(mdata.jqnumbers.shape[0], 0)
        mdata.jqnumbers = np.array([], dtype=np.float64)

        mdata.jvalues = 1, 12, 15, 3, 12, 4
        self.assertCountEqual(mdata.jqnumbers, np.array([1, 3, 4, 12, 15]))
        mdata.jrange = 1, 4
        self.assertCountEqual(mdata.jqnumbers, np.array([1, 2, 3, 4, 12, 15]))

        mdata.referencej = 2
        self.assertEqual(mdata.jqnumbers[0], 2)
        self.assertCountEqual(mdata.jqnumbers, set(mdata.jqnumbers))

    def test_parity_property(self):

        mdata = MoleculeData()

        mdata.parities = 1,
        self.assertEqual(len(mdata.pars), 1)

        mdata.parities = 1, 0
        self.assertCountEqual(mdata.pars, [0, 1])

        mdata.parities = 0
        self.assertEqual(len(mdata.pars), 1)
