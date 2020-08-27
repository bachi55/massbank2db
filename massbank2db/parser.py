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


def compile_regex(regex):
    return [re.compile(rx) for rx in regex]


def get_meta_regex(compile=True) -> dict:
    """ Create a dictionary of regex for extracting the meta data for the spectra
    """
    regex = {'accession': ['^ACCESSION:(.*)$'],
             'copyright': ['^COPYRIGHT:\s+(.*)'],
             'origin': ['^origin(?:=|:)(.*)$'],
             'record_title': ['^RECORD_TITLE:\s+(.*)$']}

    if compile:
        regex = {k: compile_regex(v) for k, v in regex.items()}

    return regex


def get_ms_regex(compile=True) -> dict:
    """ Create a dictionary of regex for extracting the Mass-spectra (MS) information for the spectra
    """
    regex = {'precursor_mz': ['^MS\$FOCUSED_ION:\s+PRECURSOR_M/Z\s+(\d*[.,]?\d*)$'],
             'precursor_type': ['^MS\$FOCUSED_ION:\s+PRECURSOR_TYPE\s+(.*)$'],
             'base_peak': ['^MS\$FOCUSED_ION:\s+BASE_PEAK\s+(\d*[.,]?\d*)$']}

    if compile:
        regex = {k: compile_regex(v) for k, v in regex.items()}

    return regex


def get_CH_regex(compile=True) -> dict:
    regex = {
        "name": ['^CH\$NAME:\s+(.*)$'],
        "molecular_formula": ['^CH\$FORMULA:\s+(.*)$'],
        "molecular_weight": ['^CH\$MOLECULAR_WEIGHT:\s+(.*)$'],
        "pubchem_id": ['^CH\$LINK:\s+PUBCHEM\s+CID:(.*)$'],  # TODO: Handle more formats with additional expressions.
        "exact_mass": ['^CH\$EXACT_MASS:\s+(.*)$'],
        "smiles": ['^CH\$SMILES:\s+(.*)$'],
        "inchikey": ['^CH\$LINK:\s+INCHIKEY\s+(.*)$'],
    }

    if compile:
        regex = {k: compile_regex(v) for k, v in regex.items()}

    return regex


def get_AC_regex(compile=True) -> dict:
    regex = {'instrument_type':     ['^AC\$INSTRUMENT_TYPE:\s+(.*)$'],
             'instrument':          ['^AC\$INSTRUMENT:\s+(.*)$'],
             # Mass Spectrometry
             'collision_energy':    ['^AC\$MASS_SPECTROMETRY:\s+COLLISION_ENERGY\s+(.*)$'],
             'ms_type':             ['^AC\$MASS_SPECTROMETRY:\s+MS_TYPE\s+(MS\d*)$'],
             'resolution':          ['^AC\$MASS_SPECTROMETRY:\s+RESOLUTION\s+(.*)$'],
             'ion_mode':            ['^AC\$MASS_SPECTROMETRY:\s+ION_MODE\s+(.*)$'],
             'fragmentation_type':  ['^AC\$MASS_SPECTROMETRY:\s+FRAGMENTATION_MODE\s+(.*)$'],
             'mass_accuracy':       ['^AC\$MASS_SPECTROMETRY:\s+ACCURACY\s+(.*)$'],
             'mass_error':          ['^AC\$MASS_SPECTROMETRY:\s+ERROR\s+(.*)$'],
             # Chromatography
             'column_name':         ['^AC\$CHROMATOGRAPHY:\s+COLUMN.*NAME\s(.*)$'],
             'flow_gradient':       ['^AC\$CHROMATOGRAPHY:\s+FLOW.*GRADIENT\s(.*)$'],
             'flow_rate':           ['^AC\$CHROMATOGRAPHY:\s+FLOW.*RATE\s(.*)$'],
             'retention_time':      ['^AC\$CHROMATOGRAPHY:\s+RETENTION.*TIME\s+(\d+[.,]?\d*)\s*($|min|sec)'],
             'solvent_A':           ['^AC\$CHROMATOGRAPHY:\s+SOLVENT\sA\s(.*)$'],
             'solvent_B':           ['^AC\$CHROMATOGRAPHY:\s+SOLVENT\sB\s(.*)$'],
             'solvent':             ['^AC\$CHROMATOGRAPHY:\s+SOLVENT\s(.*)$'],
             'column_temperature':  ['^AC\$CHROMATOGRAPHY:\s+COLUMN.*TEMPERATURE\s(.*)$']}

    if compile:
        regex = {k: compile_regex(v) for k, v in regex.items()}

    return regex


def parse_info(msfn: str, regex: dict) -> (dict, dict):
    """Parse and extract all meta data by looping through the dictionary of meta_info regexs

    Args:
         line (str): line of the msp file
    """
    infos = {key: None for key in regex}

    with open(msfn, "r") as msfile:
        line = msfile.readline()
        while not line.startswith("PK$NUM_PEAK:"):
            if _extract_info_from_line(line, regex, infos):
                continue

            line = msfile.readline()

    return infos


def _extract_info_from_line(line: str, regex: dict, out: dict) -> bool:
    for info in regex:
        if not out[info]:
            match = re.search(regex[info], line)  # re.IGNORECASE
            if match:
                out[info] = tuple(map(str.strip, match.groups()))
                return True

    return False


def parse_peaks(msfn: str):
    """

    :param msfn:
    :return:
    """
    # TODO: Check whether we are parsing peaks or annotations

    with open(msfn, "r") as msfile:
        line = msfile.readline()
        while not line.startswith("PK$NUM_PEAK:"):
            line = msfile.readline()

        # Extract peaks
        num_peaks = int(line[len("PK$NUM_PEAK: "):].strip())
        peak_list = []
        _ = msfile.readline()  # skip: PK$PEAK: m/z int. rel.int.
        line = msfile.readline()
        while line != "//\n":
            match = re.match("(\d*[,.]?\d*) (\d*[,.]?\d*) (\d*)", line.strip())
            assert match
            peak_list.append(match.groups())

            line = msfile.readline()
        assert len(peak_list) == num_peaks, "Length of extracted peak list must be equal 'NUM_PEAK'."

    return peak_list


if __name__ == "__main__":
    msfn = "/run/media/bach/EVO500GB/data/MassBank/Chubu_Univ/UT001973.txt"

    # Read all meta information from the MS-file
    info = parse_info(msfn, regex={**get_meta_regex(), **get_AC_regex()})
    print(info)

    # Read peaklist
    peaks = parse_peaks(msfn)
    print(peaks)