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

from typing import Dict, List

from massbank2db.parser import get_meta_regex, get_AC_regex, get_CH_regex, get_ms_regex


class MBSpectrum(object):
    def __init__(self, fn):
        self._fn = fn
        self._mb_raw = self._read_mb_file(self._fn)
        self._meta_information, i = self._parse_meta_information(
            self._mb_raw, {**get_meta_regex(), **get_ms_regex(), **get_CH_regex(), **get_AC_regex()})
        self._peaks = self._parse_peaks(self._mb_raw[i:])

    @staticmethod
    def _read_mb_file(fn):
        with open(fn, "r") as mbfile:
            lines = mbfile.readlines()
        return lines

    @staticmethod
    def _parse_meta_information(lines, regex: Dict[str, List[re.Pattern]]):
        meta = {k: None for k in regex}
        i = 0
        while not lines[i].startswith("PK$NUM_PEAK:"):
            for information, rxs in regex.items():
                if not meta[information]:
                    # Match the different patterns for the given information
                    match = None
                    j = 0
                    while (not match) and (j < len(rxs)):
                        match = rxs[j].match(lines[i])
                        j += 1

                    # If we have found the information, we can stop processing this line
                    if match:
                        meta[information] = tuple(map(str.strip, match.groups()))
                        break
            i += 1

        return meta, i

    @staticmethod
    def _parse_peaks(lines):
        assert lines[0].startswith("PK$NUM_PEAK:")

        num_peaks = int(lines[0][len("PK$NUM_PEAK: "):].strip())

        # Extract peaks
        peak_list = []
        i = 2  # skip: PK$PEAK: m/z int. rel.int.
        peak_regex = re.compile("(\d*[,.]?\d*) (\d*[,.]?\d*) (\d*)")
        while lines[i] != "//\n":
            match = peak_regex.match(lines[i].strip())
            assert match
            peak_list.append(match.groups())
            i += 1

        assert len(peak_list) == num_peaks, "Length of extracted peak list must be equal 'NUM_PEAK'."

        return peak_list


if __name__ == "__main__":
    msfn = "/run/media/bach/EVO500GB/data/MassBank/Chubu_Univ/UT001973.txt"

    # Read all meta information from the MS-file
    spec = MBSpectrum(msfn)

    print("hey")