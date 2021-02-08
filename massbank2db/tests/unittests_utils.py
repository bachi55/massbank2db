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

from massbank2db.utils import estimate_column_deadtime


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
