####
#
# The 'massbank2db' package can be used to an build SQLite DB from the Massbank MS/MS repository.
#
#     Copyright (C) 2020 - 2021  Eric Bach <eric.bach@aalto.fi>
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
####
import numpy as np

from typing import Optional


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


def get_mass_from_precursor_mz(precursor_mz, precursor_type):
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


def get_mass_error_in_ppm(mass, exp_precursor_mz, precursor_type):
    """
    :param mass: scalar, mass of a compound, e.g. monoisotopic or exact mass.

    :param exp_precursor_mz: scalar, ion mass / precursor mz

    :param precursor_type: string, precursor type, e.g. '[M+H]+'

    :return: scalar, mass error in ppm
    """
    try:
        calc_precursor_mz = get_precursor_mz(mass, precursor_type)
    except KeyError:
        # Precursor type was not found
        return None

    return _get_mass_error_in_ppm(calc_precursor_mz, exp_precursor_mz)


def _get_mass_error_in_ppm(calc_mass, exp_mass):
    """

    :param calc_mass:
    :param exp_mass:
    :return:
    """
    return (exp_mass - calc_mass) / calc_mass * 1e6


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

