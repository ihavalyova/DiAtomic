from diatom import *
import unittest

class TestMass(unittest.TestCase):

    def test_hydrogen_reduced_mass_value(self):

        mdata = MoleculeData()
        self.assertLess(
            mdata._calculate_molecule_reduced_mass('1H1H'), 1
        )

    def test_molecule_symbol_fail(self):

        mdata = MoleculeData()
        self.assert_raised_error(mdata, 'H1H1')
        self.assert_raised_error(mdata, '1HH')
        self.assert_raised_error(mdata, 'H1H1H')

    def assert_raised_error(self, mdata, symbol):

        self.assertRaises(
            SystemExit, mdata._calculate_molecule_reduced_mass, symbol
        )
    
    def test_molecule_symbol_fail_integer(self):

        mdata = MoleculeData()
        self.assertRaises(
            AttributeError, mdata._calculate_molecule_reduced_mass, 1
        )
    
    def test_molecule_symbol_pass(self):

        mdata = MoleculeData()
        self.assert_is_instance_of_float(mdata, '1H 1H')
        self.assert_is_instance_of_float(mdata, '1 H 1 H')
        self.assert_is_instance_of_float(mdata, '1 H 1H')
        self.assert_is_instance_of_float(mdata, '1H 1 H')
        self.assert_is_instance_of_float(mdata, ' 1H1H')
        self.assert_is_instance_of_float(mdata, '1H1H ')
    
    def assert_is_instance_of_float(self, mdata, symbol):
        self.assertIsInstance(
            mdata._calculate_molecule_reduced_mass(symbol), float
        )
