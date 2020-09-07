####
#
# The MIT License (MIT)
#
# Copyright 2020 Eric Bach <eric.bach@aalto.fi>
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
import unittest

from massbank2db.spectrum import MBSpectrum


class TestMBSpectrum(unittest.TestCase):
    def test_peak_parsing(self):
        # Spectrum 1
        peaks = [(65.0388, 454422.7),
                 (105.0703, 3594629.8),
                 (115.0544, 916283.5),
                 (121.0287, 482029.8),
                 (139.0544, 904286),
                 (163.0549, 764499.2),
                 (164.0627, 2677649.4),
                 (165.0706, 293933491.2),
                 (166.0782, 235412982.7),
                 (183.0812, 3256920.2),
                 (193.0767, 1661789.6),
                 (199.0315, 5538869.4),
                 (201.0474, 7384944.1)]
        spec = MBSpectrum("./EQ308406.txt")
        self.assertEqual(peaks, spec.get_peak_list_as_tuples())

        # Spectrum 2
        peaks = [
            (116.050500, 2998.000000),
            (117.054700, 236.000000),
            (118.029100, 170.000000),
            (131.037600, 1.23552e13),
            (132.045500, 1241.000000),
            (133.048300, 116.000000),
            (159.032400, 1732.000000),
            (160.040200, 10392.000000),
            (161.043200, 1399.000000)]
        spec = MBSpectrum("./FIO00665.txt")
        self.assertEqual(peaks, spec.get_peak_list_as_tuples())


if __name__ == '__main__':
    unittest.main()
