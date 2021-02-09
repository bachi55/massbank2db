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
import numpy as np

from typing import Optional

# TODO: Currently unsupported precursor-types
#  '[M+CH3COOH-H]-': None,
#  '[M-H-C6H10O5]-': None,
#  '[M-H+CH2O2]-': None,
#  '[M-H+C2H2O]-': None,
#  '[M-CH3]-': None,
#  '[M+H-NH3]+': None,
#  '[M+H-C9H10O5]+': None,
#  '[M-C6H10O5+H]+': None

# Values are taken from:
# https://docs.google.com/spreadsheets/d/1r4dPw1shIEy_W2BkfgPsihinwg-Nah654VlNTn8Gxo0/edit?usp=sharing
ADDUCT_MASSES = {
    # POSITIVE
    '[M+H]+': 1.007276,
    '[M-H2O+H]+': 1.007276 - 18.010565, '[M+H-H2O]+': 1.007276 - 18.010565,   # exact mass of H2O = 18.010565
    '[M-2H2O+H]+': 1.007276 - 2 * 18.010565, '[M+H-2H2O]+': 1.007276 - 2 * 18.010565,
    '[M+Na]+': 22.98922,
    '[M]+': 0,
    '[M+NH4]+': 18.03383,
    '[M+H-NH3]+': 1.007276 - 17.026549,
    # NEGATIVE
    '[M-H]-': -1.007276,
    '[M+CH3COO]-': 59.01385, '[M+CH3COOH-H]-': 59.01385,
    '[M+CHO2]-': 44.9982, '[M+HCOO]-': 44.9982, '[M+CH2O2-H]-': 44.9982, '[M-H+CH2O2]-': 44.9982,
    '[M]-': 0,
    '[M-2H]-': -1.007276
}


def get_precursor_mz(mass, precursor_type):
    """
    Calculate precursor mz based on (exact or monoisotopic) mass and precursor type

    :param mass: scalar, mass of a compound, e.g. monoisotopic or exact mass.

    :param precursor_type: string, precursor type, e.g. '[M+H]+'

    :return: scalar, ion mass / precursor mz
    """
    try:
        return mass + ADDUCT_MASSES[precursor_type]
    except KeyError:
        raise KeyError("Unsupported precursor-type '%s'." % precursor_type)


def get_mass_from_ion(precursor_mz, precursor_type):
    """
    Calculate monoisotopic mass from a given ion mass / precursor mz based on the precursor type.

    :param precursor_mz: scalar, ion mass / precursor mz

    :param precursor_type: string, precursor type, e.g. '[M+H]+'

    :return: scalar, monoisotopic mass of the compound
    """
    try:
        return precursor_mz - ADDUCT_MASSES[precursor_type]
    except KeyError:
        raise KeyError("Unsupported precursor-type '%s'." % precursor_type)


def get_mass_error_in_ppm(monoisotopic_mass, precursor_mz, precursor_type):
    """

    :param monoisotopic_mass:
    :param precursor_mz:
    :param precursor_type:
    :return:
    """
    try:
        theoretical_precursor_mz = get_precursor_mz(monoisotopic_mass, precursor_type)
    except KeyError:
        return None

    return (abs(theoretical_precursor_mz - precursor_mz) * 1e6) / theoretical_precursor_mz


def estimate_column_deadtime(length, diameter, flow_rate, flow_rate_unit="uL/min") -> Optional[float]:
    """
    Estimates the column's deadtime given it's length, diameter and flowrate:

        dt = ((pi * (diameter * 0.5)^2 * length * 0.5) / 1000) / flowrate(in mL/min)

    :param length: scalar, length of the column in mm
    :param diameter: scalar, diameter of the column in mm
    :param flow_rate: scalar, flowrate of the column in 'flowrate_unit'
    :param flow_rate_unit: string, defining the unit of the given flowrate
    :return: scalar, estimated column deadtime or None if for any of the required input parameters np.isnan is True.
    """
    if length is None or diameter is None or flow_rate is None:
        return None

    if flow_rate_unit.lower() == "ul/min":
        flow_rate /= 1000.0  # convert to mL/min
    elif flow_rate_unit.lower() == "ml/min":
        pass
    else:
        raise ValueError("Invalid flowrate unit: '%s'." % flow_rate_unit)

    return ((np.pi * (diameter / 2)**2 * length * 0.5) / 1000) / flow_rate


def named_row_factory(cursor, row):
    """
    SOURCE: https://docs.python.org/3.5/library/sqlite3.html#sqlite3.Connection.row_factory
    """
    return {col[0]: row[idx] for idx, col in enumerate(cursor.description)}


def get_ppm_window(mass, ppm):
    abs_deviation = mass / 1e6 * ppm
    return mass - abs_deviation, mass + abs_deviation

