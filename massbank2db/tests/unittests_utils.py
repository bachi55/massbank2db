####
#
# The MIT License (MIT)
#
# Copyright 2021 Eric Bach <eric.bach@aalto.fi>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
####
import numpy as np
import unittest

from massbank2db.utils import estimate_column_deadtime, get_precursor_mz, get_mass_from_precursor_mz, ADDUCT_MASSES
from massbank2db.utils import _get_mass_error_in_ppm, get_ppm_window


class TestMassTools(unittest.TestCase):
    def test_calculate_precursor_mz(self):
        # Test against examples from: https://onlinelibrary.wiley.com/doi/full/10.1002/rcm.3102
        # Atrazine
        self.assertEqual(216.1010, np.round(get_precursor_mz(mass=215.093773, precursor_type="[M+H]+"), 4))
        # Carbofuran
        self.assertEqual(222.1125, np.round(get_precursor_mz(mass=221.105193, precursor_type="[M+H]+"), 4))
        # Diphenhydramine [M+H]+
        self.assertEqual(256.1696, np.round(get_precursor_mz(mass=255.162314, precursor_type="[M+H]+"), 4))
        # Diphenhydramine [M+Na]+
        self.assertEqual(255.1623 + 22.9892, np.round(get_precursor_mz(mass=255.162314, precursor_type="[M+Na]+"), 4))
        # Naproxen
        self.assertEqual(229.0870, np.round(get_precursor_mz(mass=230.094294, precursor_type="[M-H]-"), 4))

    def test_conversion_from_mass_to_ion_and_back(self):
        precursor_types = list(ADDUCT_MASSES.keys())

        for rep in range(100):
            mass = np.random.RandomState(rep).rand() * 400 + 100
            precursor_type = np.random.RandomState(rep).choice(precursor_types)
            self.assertAlmostEqual(mass,
                                   get_mass_from_precursor_mz(get_precursor_mz(mass, precursor_type), precursor_type))

    def test_mass_error_in_ppm(self):
        # Test against examples from: https://onlinelibrary.wiley.com/doi/full/10.1002/rcm.3102
        self.assertEqual(2.8, np.round(np.abs(_get_mass_error_in_ppm(216.1016, 216.1010)), 1))
        self.assertEqual(0, _get_mass_error_in_ppm(216.1010, 216.1010))
        self.assertEqual(2.0, np.round(np.abs(_get_mass_error_in_ppm(250.0190, 250.0185)), 1))

        # Test against: https://warwick.ac.uk/fac/sci/chemistry/research/barrow/barrowgroup/calculators/mass_errors/
        self.assertEqual(np.round(1.726473, 4), np.abs(np.round(_get_mass_error_in_ppm(463.37243, 463.37323), 4)))
        self.assertEqual(np.round(-1221.001221, 4), np.round(_get_mass_error_in_ppm(819, 818), 4))
        self.assertEqual(np.round(-12.208373, 4), np.round(_get_mass_error_in_ppm(819.11, 819.1), 4))

    def test_ppm_window(self):
        # Test against: https://warwick.ac.uk/fac/sci/chemistry/research/barrow/barrowgroup/calculators/mass_errors/
        self.assertEqual((399.9996, 400.0004), get_ppm_window(400, 1))
        self.assertEqual((399.998, 400.002), get_ppm_window(400, 5))
        _min, _max = get_ppm_window(271.4324, 5)
        self.assertEqual((271.431043, 271.433757), (np.round(_min, 6), np.round(_max, 6)))


class TestColumnDeadtimeEstimation(unittest.TestCase):
    def test_bordercases(self):
        self.assertEqual(estimate_column_deadtime(None, 1, 1), None)
        self.assertEqual(estimate_column_deadtime(1, None, 1), None)
        self.assertEqual(estimate_column_deadtime(1, 1, None), None)

    def test_correctness(self):
        np.testing.assert_allclose(estimate_column_deadtime(150, 2.1, 0.5, flow_rate_unit="mL/Min"),
                                   0.519540885087412, )
        np.testing.assert_allclose(estimate_column_deadtime(150, 2.1, 500, flow_rate_unit="uL/Min"),
                                   0.519540885087412)
        np.testing.assert_allclose(estimate_column_deadtime(150, 2.1, 0.25, flow_rate_unit="mL/Min"),
                                   1.03908177017482)
        np.testing.assert_allclose(estimate_column_deadtime(150, 2.1, 250, flow_rate_unit="uL/Min"),
                                   1.03908177017482)
        np.testing.assert_allclose(estimate_column_deadtime(100, 2.1, 0.3, flow_rate_unit="mL/Min"),
                                   0.577267650097125)
        np.testing.assert_allclose(estimate_column_deadtime(100, 2.1, 300, flow_rate_unit="uL/Min"),
                                   0.577267650097125)


if __name__ == '__main__':
    unittest.main()
