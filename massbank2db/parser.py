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


def _compile_regex(regex):
    return [re.compile(rx) for rx in regex]


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
    regex = {'accession':    ['^ACCESSION:\s+(.*)$'],
             'deprecated':   ['^DEPRECATED:\s+(.*)$'],
             'copyright':    ['^COPYRIGHT:\s+(.*)'],
             'origin':       ['^origin(?:=|:)(.*)$'],
             'record_title': ['^RECORD_TITLE:\s+(.*)$'],
             'license':      ['^LICENSE:\s+(.*)$']}

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
    regex = {'precursor_mz':    ['^MS\$FOCUSED_ION:\s+PRECURSOR_M/Z\s+(\d*[.,]?\d*)$'],
             'precursor_type':  ['^MS\$FOCUSED_ION:\s+PRECURSOR_TYPE\s+(.*)$'],
             'base_peak':       ['^MS\$FOCUSED_ION:\s+BASE_PEAK\s+(\d*[.,]?\d*)$']}

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
        "name":              ['^CH\$NAME:\s+(.*)$'],
        "molecular_formula": ['^CH\$FORMULA:\s+(.*)$'],
        "molecular_weight":  ['^CH\$MOLECULAR_WEIGHT:\s+(.*)$'],
        # As PubChem ID we consider only the CIDs (not SIDs!)
        "pubchem_id":        ['^CH\$LINK:\s+PUBCHEM\s+CID:(\d+)$',
                              '^CH\$LINK:\s+PUBCHEM\s+CID:(\d+)\s+SID:\d+$',
                              '^CH\$LINK:\s+PUBCHEM\s+SID:\d+\s+CID:(\d+)$'],
        "exact_mass":        ['^CH\$EXACT_MASS:\s+(.*)$'],
        # MassBank documentation states, that the SMILES strings correspond to isomeric smiles:
        # https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#225-chsmiles
        "smiles_iso":        ['^CH\$SMILES:\s+(.*)$'],
        "inchikey":          ['^CH\$LINK:\s+INCHIKEY\s+(.*)$'],
        "inchi":             ['^CH\$IUPAC:\s+(.*)$']
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
    regex = {'instrument_type':     ['^AC\$INSTRUMENT_TYPE:\s+(.*)$'],
             'instrument':          ['^AC\$INSTRUMENT:\s+(.*)$'],
             # Mass Spectrometry
             'collision_energy':    ['^AC\$MASS_SPECTROMETRY:\s+COLLISION_ENERGY\s+(.*)$'],
             'ms_type':             ['^AC\$MASS_SPECTROMETRY:\s+MS_TYPE\s+(MS\d*)$'],
             'resolution':          ['^AC\$MASS_SPECTROMETRY:\s+RESOLUTION\s+(.*)$'],
             'ion_mode':            ['^AC\$MASS_SPECTROMETRY:\s+ION_MODE\s+(.*)$'],
             'fragmentation_mode':  ['^AC\$MASS_SPECTROMETRY:\s+FRAGMENTATION_MODE\s+(.*)$'],
             'mass_accuracy':       ['^AC\$MASS_SPECTROMETRY:\s+ACCURACY\s+(.*)$'],
             'mass_error':          ['^AC\$MASS_SPECTROMETRY:\s+ERROR\s+(.*)$'],
             # Chromatography
             'column_name':         ['^AC\$CHROMATOGRAPHY:\s+COLUMN.*NAME\s(.*)$'],
             'flow_gradient':       ['^AC\$CHROMATOGRAPHY:\s+FLOW.*GRADIENT\s(.*)$'],
             'flow_rate':           ['^AC\$CHROMATOGRAPHY:\s+FLOW.*RATE\s(.*)$'],
             'retention_time':      ['^AC\$CHROMATOGRAPHY:\s+RETENTION.*TIME\s+(\d+[.,]?\d*)\s*($|min|sec|s|m)'],
             'solvent_A':           ['^AC\$CHROMATOGRAPHY:\s+SOLVENT\sA\s(.*)$'],
             'solvent_B':           ['^AC\$CHROMATOGRAPHY:\s+SOLVENT\sB\s(.*)$'],
             'solvent':             ['^AC\$CHROMATOGRAPHY:\s+SOLVENT\s(.*)$'],
             'column_temperature':  ['^AC\$CHROMATOGRAPHY:\s+COLUMN.*TEMPERATURE\s(.*)$']}

    if compile:
        regex = {k: _compile_regex(v) for k, v in regex.items()}

    return regex


if __name__ == "__main__":
    import os

    from glob import glob

    from massbank2db.spectrum import MBSpectrum

    # msfn = "/run/media/bach/EVO500GB/data/MassBank/Chubu_Univ/UT001973.txt"

    for msfn in sorted(glob("/run/media/bach/EVO500GB/data/MassBank-data_bachi55/ISAS_Dortmund/IA*.txt")):
        # Read all meta information from the MS-file
        try:
            spec = MBSpectrum(msfn)
        except AssertionError:
            print(os.path.basename(msfn))
