####
#
# The MIT License (MIT)
#
# Copyright 2020, 2021 Eric Bach <eric.bach@aalto.fi>
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

# TODO: Currently unsupported precursor-types
#  '[M+CH3COOH-H]-': None,
#  '[M-2H2O+H]+': None,
#  '[M-H-C6H10O5]-': None,
#  '[M-H+CH2O2]-': None,
#  '[M-H+C2H2O]-': None,
#  '[M-CH3]-': None,
#  '[M-H2O+H]+': None,
#  '[M+H-H2O]+': None,
#  '[M+H-NH3]+': None,
#  '[M+H-C9H10O5]+': None,
#  '[M-C6H10O5+H]+': None

def get_precursor_mz(exact_mass, precursor_type):
    """
    Calculate precursor mz based on exact mass and precursor type

    :param exact_mass: float, exact mass of compound of interest

    :param precursor_type: string, precursor type, e.g. [M+H]+

    :return: scalar, neutral mass of compound
    """
    # Values are taken from:
    # https://docs.google.com/spreadsheets/d/1r4dPw1shIEy_W2BkfgPsihinwg-Nah654VlNTn8Gxo0/edit?usp=sharing
    d = {'[M+H]+': 1.007276,
         '[M-H]-': -1.007276,
         '[M+HCOO]-': -44.9982,
         '[M+CH3COO]-': -59.01385,
         '[M+Na]+': 22.98922,
         '[M]+': 0,
         '[M]-': 0,
         '[M-2H]-': -1.00728,
         '[M+NH4]+': 18.03383}

    try:
        return exact_mass + d[precursor_type]
    except KeyError:
        raise KeyError("Unsupported precursor-type '%s'." % precursor_type)


def get_mass_error_in_ppm(exact_mass, precursor_mz, precursor_type):
    """

    :param exact_mass:
    :param precursor_mz:
    :param precursor_type:
    :return:
    """
    try:
        theoretical_precursor_mz = get_precursor_mz(exact_mass, precursor_type)
    except KeyError:
        return None

    return (abs(theoretical_precursor_mz - precursor_mz) * 1e6) / theoretical_precursor_mz
