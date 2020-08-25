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


def get_meta_regex() -> dict:
    """ Create a dictionary of regex for extracting the meta data for the spectra
    """
    # NOTE: will just ignore cases, to avoid repetition here
    meta_parse = {'collision_energy': '^AC\$MASS_SPECTROMETRY:\s+COLLISION_ENERGY\s+(.*)$',
                  'ms_level': '^AC\$MASS_SPECTROMETRY:\s+MS_TYPE\s+\D*(\d*)$',
                  'accession': '^ACCESSION:(.*)$',
                  'resolution': '^AC\$MASS_SPECTROMETRY:\s+RESOLUTION\s+(.*)$',
                  'ion_mode': '^AC\$MASS_SPECTROMETRY:\s+ION_MODE\s+(.*)$',
                  'fragmentation_type': '^AC\$MASS_SPECTROMETRY:\s+FRAGMENTATION_MODE\s+(.*)$',
                  'instrument_type': '^AC\$INSTRUMENT_TYPE:\s+(.*)$',
                  'instrument': '^AC\$INSTRUMENT:\s+(.*)$',
                  'copyright': '^COPYRIGHT:\s+(.*)',
                  'mass_accuracy': '^AC\$MASS_SPECTROMETRY:\s+ACCURACY\s+(.*)$',
                  'mass_error': '^AC\$MASS_SPECTROMETRY:\s+ERROR\s+(.*)$',
                  'origin': '^origin(?:=|:)(.*)$',
                  'record_title': '^RECORD_TITLE:\s+(.*)$'}

    return meta_parse


def get_ms_regex() -> dict:
    """ Create a dictionary of regex for extracting the Mass-spectra (MS) information for the spectra
    """
    meta_parse = {'precursor_mz': '^MS\$FOCUSED_ION:\s+PRECURSOR_M/Z\s+(\d*[.,]?\d*)$',
                  'precursor_type': '^MS\$FOCUSED_ION:\s+PRECURSOR_TYPE\s+(.*)$'}

    return meta_parse


def get_chromatographic_regex() -> dict:
    """ Create a dictionary of regex for extracting the Chromatographic (AC$CHROMATOGRAPHY) information for the spectra
    """
    meta_parse = {
        'column_name': '^AC\$CHROMATOGRAPHY:\s+COLUMN.*NAME\s(.*)$',
        'flow_gradient': '^AC\$CHROMATOGRAPHY:\s+FLOW.*GRADIENT\s(.*)$',
        'flow_rate': '^AC\$CHROMATOGRAPHY:\s+FLOW.*RATE\s(.*)$',
        'retention_time': '^AC\$CHROMATOGRAPHY:\s+RETENTION.*TIME\s+(\d+[.,]?\d*)\s?($|min|sec)',
        'solvent_A': '^AC\$CHROMATOGRAPHY:\s+SOLVENT\sA\s(.*)$',
        'solvent_B': '^AC\$CHROMATOGRAPHY:\s+SOLVENT\sB\s(.*)$',
        'solvent': '^AC\$CHROMATOGRAPHY:\s+SOLVENT\s(.*)$',
        'column_temperature': '^AC\$CHROMATOGRAPHY:\s+COLUMN.*TEMPERATURE\s(.*)$'
    }

    return meta_parse


def get_compound_regex() -> dict:
    """ Create a dictionary of regex for extracting the compound information for the spectra
    """

    # NOTE: will just ignore cases in the regex, to avoid repetition here
    meta_parse = {'name': '^CH\$NAME:\s+(.*)$',
                  'other_names': '^CH\$NAME:\s+(.*)$',
                  'inchikey_id': '^CH\$LINK:\s+INCHIKEY\s+(.*)$',
                  'molecular_formula': '^CH\$FORMULA:\s+(.*)$',
                  'molecular_weight': '^CH\$MOLECULAR_WEIGHT:\s+(.*)$',
                  'pubchem_id': '^CH\$LINK:\s+PUBCHEM\s+CID:(.*)$',
                  'chemspider_id': '^CH\$LINK:\s+CHEMSPIDER\s+(.*)$',
                  'compound_class': '^CH\$COMPOUND_CLASS:\s+(.*)$',
                  'exact_mass': '^CH\$EXACT_MASS:\s+(.*)$',
                  'smiles': '^CH\$SMILES:\s+(.*)$'}

    return meta_parse


def parse_info(msfn: str) -> (dict, dict):
    """Parse and extract all meta data by looping through the dictionary of meta_info regexs

    updates self.meta_info

    Args:
         line (str): line of the msp file
    """
    # Get meta information patterns to search for
    meta_regex = get_meta_regex()
    meta_infos = {key: None for key in meta_regex}

    # Get MS information patterns to search for
    ms_regex = get_ms_regex()
    ms_infos = {key: None for key in ms_regex}

    # Get chromatographic information patterns to search for
    chromatographic_regex = get_chromatographic_regex()
    chromatographic_infos = {key: None for key in chromatographic_regex}

    # Get compound information patterns to search for
    cmp_regex = get_compound_regex()
    cmp_infos = {key: None for key in cmp_regex}

    with open(msfn, "r") as msfile:
        line = msfile.readline()
        while not line.startswith("PK$NUM_PEAK:"):
            if _extract_info_from_line(line, meta_regex, meta_infos):
                continue

            if _extract_info_from_line(line, ms_regex, ms_infos):
                continue

            if _extract_info_from_line(line, chromatographic_regex, chromatographic_infos):
                continue

            if _extract_info_from_line(line, cmp_regex, cmp_infos):
                continue

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

    return meta_infos, ms_infos, cmp_infos, chromatographic_infos, peak_list


def _extract_info_from_line(line: str, regex: dict, out: dict) -> bool:
    for info in regex:
        if not out[info]:
            match = re.search(regex[info], line, re.IGNORECASE)
            if match:
                out[info] = tuple(map(str.strip, match.groups()))

                return True


if __name__ == "__main__":
    msfn = "/run/media/bach/EVO500GB/data/MassBank/Chubu_Univ/UT001973.txt"

    # Read all lines in the MS-file
    for info in parse_info(msfn):
        print(info)
