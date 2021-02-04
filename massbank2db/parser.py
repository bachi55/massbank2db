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
import re

from typing import Union, List, Optional, Callable, Tuple


def _compile_regex(regex: Union[str, List[str]], *args) -> Union[re.Pattern, List[re.Pattern]]:
    if isinstance(regex, str):
        return re.compile(regex, *args)
    elif isinstance(regex, list):
        return [re.compile(rx, *args) for rx in regex]
    else:
        raise ValueError("Invalid input. Regular expressions need to be provided as list[str] or str.")


def get_meta_regex(compile=True) -> dict:
    """
    Create a dictionary of regex for extracting the meta data for the spectra

    See also: https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#21-record-specific-information

    :param compile: boolean indicating whether the regular expressions should be compiled already.

    :return: dictionary,
        keys, strings represent to name of the property to extract (e.g. same name as in the MBSpectrum class)
        values,
            list of strings, regular expressions to extract the particular information

            or (if compile=True)

            list of re.Pattern, pre-compiled regular expressions to extract the information
    """
    regex = {'accession':    [r'^ACCESSION:\s+(.*)$'],
             'deprecated':   [r'^DEPRECATED:\s+(.*)$'],
             'copyright':    [r'^COPYRIGHT:\s+(.*)'],
             'origin':       [r'^origin(?:=|:)(.*)$'],
             'record_title': [r'^RECORD_TITLE:\s+(.*)$'],
             'license':      [r'^LICENSE:\s+(.*)$']}

    if compile:
        regex = {k: _compile_regex(v) for k, v in regex.items()}

    return regex


def get_ms_regex(compile=True) -> dict:
    """
    Create a dictionary of regex for extracting the Mass-spectra (MS) information for the spectra.

    See also: https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#25-description-of-mass-spectral-data

    :param compile: boolean indicating whether the regular expressions should be compiled already.

    :return: dictionary,
        keys, strings represent to name of the property to extract (e.g. same name as in the MBSpectrum class)
        values,
            list of strings, regular expressions to extract the particular information

            or (if compile=True)

            list of re.Pattern, pre-compiled regular expressions to extract the information
    """
    regex = {'precursor_mz':    [r'^MS\$FOCUSED_ION:\s+PRECURSOR_M/Z\s+(\d*[.,]?\d*)$'],
             'precursor_type':  [r'^MS\$FOCUSED_ION:\s+PRECURSOR_TYPE\s+(.*)$'],
             'base_peak':       [r'^MS\$FOCUSED_ION:\s+BASE_PEAK\s+(\d*[.,]?\d*)$']}

    if compile:
        regex = {k: _compile_regex(v) for k, v in regex.items()}

    return regex


def get_CH_regex(compile=True) -> dict:
    """
    Regular expressions related to the ground truth molecular structure annotation of the spectra. The information are
    sufficient to represent the structure of the molecule.

    See also: https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#22-information-of-chemical-compound-analyzed

    :param compile: boolean indicating whether the regular expressions should be compiled already.

    :return: dictionary,
        keys, strings represent to name of the property to extract (e.g. same name as in the MBSpectrum class)
        values,
            list of strings, regular expressions to extract the particular information

            or (if compile=True)

            list of re.Pattern, pre-compiled regular expressions to extract the information
    """
    regex = {
        "name":              [r'^CH\$NAME:\s+(.*)$'],
        "molecular_formula": [r'^CH\$FORMULA:\s+(.*)$'],
        "molecular_weight":  [r'^CH\$MOLECULAR_WEIGHT:\s+(.*)$'],
        # As PubChem ID we consider only the CIDs (not SIDs!)
        "pubchem_id":        [r'^CH\$LINK:\s+PUBCHEM\s+CID:(\d+)$',
                              r'^CH\$LINK:\s+PUBCHEM\s+CID:(\d+)\s+SID:\d+$',
                              r'^CH\$LINK:\s+PUBCHEM\s+SID:\d+\s+CID:(\d+)$'],
        # Massbank documentation states, that the exact mass field stores the 'monoisotopic mass'
        # https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#224-chexact_mass
        "monoisotopic_mass": [r'^CH\$EXACT_MASS:\s+(.*)$'],
        # MassBank documentation states, that the SMILES strings correspond to isomeric smiles:
        # https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#225-chsmiles
        "smiles_iso":        [r'^CH\$SMILES:\s+(.*)$'],
        "inchikey":          [r'^CH\$LINK:\s+INCHIKEY\s+(.*)$'],
        "inchi":             [r'^CH\$IUPAC:\s+(.*)$']
    }

    if compile:
        regex = {k: _compile_regex(v) for k, v in regex.items()}

    return regex


def get_AC_regex(compile=True) -> dict:
    """
    Regular expressions related to the analytical method and conditions of the MassBank entry.

    See also: https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#24-analytical-method-and-conditions

    :param compile: boolean indicating whether the regular expressions should be compiled already.

    :return: dictionary,
        keys, strings represent to name of the property to extract (e.g. same name as in the MBSpectrum class)
        values,
            list of strings, regular expressions to extract the particular information

            or (if compile=True)

            list of re.Pattern, pre-compiled regular expressions to extract the information
    """
    regex = {'instrument_type':     [r'^AC\$INSTRUMENT_TYPE:\s+(.*)$'],
             'instrument':          [r'^AC\$INSTRUMENT:\s+(.*)$'],
             # Mass Spectrometry
             'collision_energy':    [r'^AC\$MASS_SPECTROMETRY:\s+COLLISION_ENERGY\s+(.*)$'],
             'ms_type':             [r'^AC\$MASS_SPECTROMETRY:\s+MS_TYPE\s+(MS\d*)$'],
             'resolution':          [r'^AC\$MASS_SPECTROMETRY:\s+RESOLUTION\s+(.*)$'],
             'ion_mode':            [r'^AC\$MASS_SPECTROMETRY:\s+ION_MODE\s+(.*)$'],
             'fragmentation_mode':  [r'^AC\$MASS_SPECTROMETRY:\s+FRAGMENTATION_MODE\s+(.*)$'],
             'mass_accuracy':       [r'^AC\$MASS_SPECTROMETRY:\s+ACCURACY\s+(.*)$'],
             'mass_error':          [r'^AC\$MASS_SPECTROMETRY:\s+ERROR\s+(.*)$'],
             # Chromatography
             'column_name':         [r'^AC\$CHROMATOGRAPHY:\s+COLUMN.*NAME\s(.*)$'],
             'flow_gradient':       [r'^AC\$CHROMATOGRAPHY:\s+FLOW.*GRADIENT\s(.*)$'],
             'flow_rate':           [r'^AC\$CHROMATOGRAPHY:\s+FLOW.*RATE\s(.*)$'],
             'retention_time':      [r'^AC\$CHROMATOGRAPHY:\s+RETENTION.*TIME\s+(\d+[.,]?\d*)\s*($|min|sec|s|m)'],
             'solvent_A':           [r'^AC\$CHROMATOGRAPHY:\s+SOLVENT\sA\s(.*)$'],
             'solvent_B':           [r'^AC\$CHROMATOGRAPHY:\s+SOLVENT\sB\s(.*)$'],
             'solvent':             [r'^AC\$CHROMATOGRAPHY:\s+SOLVENT\s(.*)$'],
             'column_temperature':  [r'^AC\$CHROMATOGRAPHY:\s+COLUMN.*TEMPERATURE\s(.*)$']}

    if compile:
        regex = {k: _compile_regex(v) for k, v in regex.items()}

    return regex


def parse_column_name(column_name: str) -> Tuple[Union[None, Tuple[float, str]], Union[None, Tuple[float, str]]]:
    """
    :param column_name: string, column name provided in massbank

    :return:
    """

    diameter = None
    length = None

    # --------------------------------------------------------
    # Extract information about the column length and diameter
    # --------------------------------------------------------
    for regex in _compile_regex([
        r"""
        (?P<diameter_value>\d.\d|\d)    # Match diameter values like '1.7' or '2'
        \s?                             # Sometime we need to skip a whitespace
        (?P<diameter_unit>|mm)          # Match diameter units like '' or 'mm'
        \s?                             # 
        [x\*]                           # Diameter x Length or vice versa ('*' and 'x' are used)
        \s?                             # 
        (?P<length_value>[\d]{2,3})     # Match length values like 50 or 120 (only two- and three-digit values)
        \s?                             #
        (?P<length_unit>mm)             # Match length units like 'mm' (needs to be given if length is last one)
        """,
        r"""
        (?P<length_value>[\d]{2,3})     # Match length values like 50 or 120 (only two- and three-digit values)
        \s?                             # 
        (?P<length_unit>|mm)            # Match length units like '' or 'mm'
        \s?                             # 
        [x\*]                           # Length x Diameter or vice versa ('*' and 'x' are used)
        \s?                             # 
        (?P<diameter_value>\d.\d|\d)    # Match diameter values like '1.7' or '2'
        \s?                             #
        (?P<diameter_unit>mm)           # Match diameter units like 'mm' (needs to be given if diameter is last one)
        """
    ], re.VERBOSE):
        match = regex.search(column_name)
        if match:
            diameter = (float(match.group("diameter_value")), match.group("diameter_unit"))
            length = (float(match.group("length_value")), match.group("length_unit"))

            break

    return diameter, length


def parse_flow_rate_string(flow_rate_str: str, aggregate_flow_rates: Optional[Callable[[List[float]], float]] = None) \
        -> Optional[Tuple[Union[List[float], float], str]]:
    """
    :param flow_rate_str:
    :param aggregate_flow_rates:
    :return:
    """
    # ---------------------------------------
    # Extract information about the flow rate
    # ---------------------------------------
    _fr_vals, _fr_units = [], []
    for regex in _compile_regex([
        r"""
        (?P<flowrate_value>\d{3})           # Matches three digit flow rates like '250' or '300'
        \s?                                 #
        (?P<flowrate_unit>u[Ll]\s?\/\s?min) # Matches flow rate units like 'uL/min' or 'ul / min'
        """,
        r"""
        (?P<flowrate_value>0.\d+)                           # Matches flow rates like '0.10', '0.25' or '0.2'
        \s?                                                 #
        (?P<flowrate_unit>m[Ll]\s?\/\s?min|m[Ll]\smin-1)    # Matches flow rate units like 'ml/ min' or 'mL min-1'
        """
    ], re.VERBOSE):
        for match in regex.finditer(flow_rate_str):
            _fr_vals.append(float(match.group("flowrate_value")))
            _fr_units.append(
                match.group("flowrate_unit")
                    .replace("min-1", "/min")       # bring things like min-1 --> /min
                    .replace(" ", "")               # remove white-spaces
            )

            if _fr_units[-1] != _fr_units[0]:
                print("All flow-rate units must be equal. Got: {} != {}.".format(_fr_units[-1], _fr_units[0]))
                return None

        if len(_fr_vals) > 0:
            break

    if aggregate_flow_rates is not None:
        _fr_vals = [aggregate_flow_rates(_fr_vals)]

    if len(_fr_vals) == 1:
        flow_rate = (_fr_vals[0], _fr_units[0])
    elif len(_fr_vals) > 1:
        flow_rate = (_fr_vals, _fr_units[0])
    else:
        flow_rate = None

    return flow_rate


def _run_parser():
    import os

    from glob import glob

    from massbank2db.spectrum import MBSpectrum

    for msfn in sorted(glob("/home/bach/Documents/doctoral/data/MassBank-data_bachi55/ISAS_Dortmund/IA*.txt")):
        # Read all meta information from the MS-file
        try:
            spec = MBSpectrum(msfn)
        except AssertionError:
            print(os.path.basename(msfn))


def _run_extraction():
    column_name = "XBridge C18 3.5um, 2.1x50mm, Waters"
    flow_rate = "200 uL/min at 0-3 min, 400 uL/min at 14 min, 480 uL/min at 16-19 min, 200 uL/min at 19.1-20 min"

    print(parse_column_name(column_name))
    print(parse_flow_rate_string(flow_rate))


if __name__ == "__main__":
    # _run_parser()
    _run_extraction()

