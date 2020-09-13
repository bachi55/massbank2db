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
import logging
import numpy as np

from ctypes import c_double, c_int, byref, cdll
from typing import Dict, List
from scipy.spatial.distance import pdist

from massbank2db.parser import get_meta_regex, get_AC_regex, get_CH_regex, get_ms_regex

# Setup the Logger
LOGGER = logging.getLogger(__name__)
CH = logging.StreamHandler()
FORMATTER = logging.Formatter('[%(levelname)s] %(name)s : %(message)s')
CH.setFormatter(FORMATTER)
LOGGER.addHandler(CH)


class MBSpectrum(object):
    def __init__(self, fn=None):
        if fn is not None:
            mb_raw = self._read_mb_file(fn)
            meta_information, i = self._parse_meta_information(
                mb_raw, {**get_meta_regex(), **get_ms_regex(), **get_CH_regex(), **get_AC_regex()})
            self._meta_information = self._sanitize_meta_information(meta_information)
            self._mz, self._int = self._parse_peaks(mb_raw[i:])
            self._mz = list(map(float, self._mz))
            self._int = list(map(float, self._int))
        else:
            self._meta_information = {}
            self._mz = []
            self._int = []

    def get(self, key, default=None):
        return self._meta_information.get(key, default)

    def set(self, key, value):
        self._meta_information[key] = value

    def get_peak_list_as_tuples(self):
        return list(zip(self._mz, self._int))

    def get_mz(self):
        return self._mz

    def get_int(self):
        return self._int

    def get_meta_information(self):
        return self._meta_information

    def set_mz(self, mz):
        self._mz = mz

    def set_int(self, ints):
        self._int = ints

    def update_molecule_structure_information_using_pubchem(self, db_conn):
        """
        Information of the molecular structure are extracted from PubChem. We use the CID (if provided) or the InChIKey
        provided in the Massbank file to find the compound in PubChem. The motivation is, that we have SMILES
        information coming from a single source and therefore being consistent.

        :param db_conn: sqlite.Connection, with a local PubChem DB stored in an SQLite file.

        :return: boolean, indicating whether an update could be performed.
        """
        # Get the ID, PubChem (cid) or InChIKey, to query information from the local PubChem DB
        if self.get("pubchem_id"):
            id, id_type = int(self.get("pubchem_id")), "cid"
        elif self.get("inchikey"):
            id, id_type = self.get("inchikey"), "InChIKey"
        else:
            LOGGER.info("(%s) : Cannot update molecule structure information. No ID defined, either CID or InChIKey."
                        % self.get("accession"))
            return False

        # Fetch information from local DB
        # TODO: Update also molecular weight. For that we need to rebuild the local PubChem DB.
        rows = db_conn.execute("SELECT cid, InChI, InChIKey, SMILES_CAN, SMILES_ISO, exact_mass, molecular_formula "
                               "    FROM compounds"
                               "    WHERE %s is ? "
                               "    ORDER BY cid ASC" % id_type, (id, )).fetchall()

        if len(rows) == 0:
            LOGGER.info("({}) Could not find any compound with {}={} in the local PubChem DB."
                        .format(self.get("accession"), id_type, id))
            return False
        elif len(rows) > 1:
            assert id_type == "InChIKey"
            LOGGER.warning("(%s) Multiple compounds (n=%d) matching the InChIKey (%s). Taking the one with the lowest "
                           "CID to update the compound information." % (self.get("accession"), len(rows), id))

        # Update the information
        for c, name in enumerate(["pubchem_id", "inchi", "inchikey", "smiles_can", "smiles_iso", "exact_mass",
                                  "molecular_formula"]):
            self.set(name, rows[0][c])

        return True

    @staticmethod
    def _sanitize_meta_information(meta_info_in):
        meta_info_out = {}
        for k, v in meta_info_in.items():
            if v is None:
                continue

            if len(v) == 0:
                LOGGER.warning("Empty information for {}={}".format(k, v))
                continue

            if len(v) == 1:
                meta_info_out[k] = v[0]
            else:
                if k == "retention_time":
                    meta_info_out["retention_time"] = float(v[0])
                    meta_info_out["retention_time_unit"] = v[1]
                else:
                    raise NotImplemented("Multiple outputs for information '%s'.")
        return meta_info_out

    @staticmethod
    def _read_mb_file(fn):
        """
        Read all lines of a MassBank entry.

        :param fn: string, path to the MassBank entry file

        :return: list of strings, all lines in the MassBank entry
        """
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
                    # Match the different patterns for the given information (stop after first match)
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
        mzs, ints = [], []
        i = 2  # skip: PK$PEAK: m/z int. rel.int.
        while lines[i] != "//\n":
            _mz, _int, _ = lines[i].strip().split(" ")
            mzs.append(_mz)
            ints.append(_int)
            i += 1

        assert len(mzs) == num_peaks, "Length of extracted peak list must be equal 'NUM_PEAK'."

        return mzs, ints

    @staticmethod
    def merge_spectra(spectra, eppm=5, eabs=0.001, rt_agg="min"):
        # TODO: Apply filter to spectra, e.g. normalize intensities, cut of high mass peaks (> precursor mz)

        mzs = []
        ints = []
        for spectrum in spectra:
            mzs += spectrum.get_mz()
            ints += spectrum.get_int()
        mzs = np.array(mzs)
        ints = np.array(ints)

        grouping = mzClust_hclust(mzs, eppm, eabs)
        print(grouping)

        mzs_out = []
        ints_out = []
        for peaks in grouping.values():
            mzs_out.append(np.mean(mzs[peaks]))
            ints_out.append(np.max(ints[peaks]))
        mzs_out = np.array(mzs_out)
        ints_out = np.array(ints_out)

        _idc_sorted = np.argsort(mzs_out)
        mzs_out = mzs_out[_idc_sorted].tolist()
        ints_out = ints_out[_idc_sorted].tolist()

        # TODO: Set properties of the output spectrum.
        spec_out = MBSpectrum()
        spec_out.set_mz(mzs_out)
        spec_out.set_int(ints_out)

        for info in spectra[0].get_meta_information():
            _info = set()
            for spectrum in spectra:
                _info.add(spectrum.get(info))

            if len(_info) == 1:
                spec_out.set(info, _info.pop())
            else:
                spec_out.set(info, [i for i in _info])

        if spec_out.get("retention_time"):
            assert isinstance(spec_out.get("retention_time_unit"), str)

            if rt_agg == "min":
                spec_out.set("retention_time", np.min(spec_out.get("retention_time")))
            elif rt_agg == "mean":
                spec_out.set("retention_time", np.mean(spec_out.get("retention_time")))
            elif rt_agg == "median":
                spec_out.set("retention_time", np.median(spec_out.get("retention_time")))
            else:
                ValueError("Invalid retention time aggregation function: '%s'. Choices are: 'min', 'mean' and "
                           "'median'.")


def mzClust_hclust(mzs, eppm, eabs):
    N = len(mzs)

    # Empty spectrum or single peak: no grouping needed
    if N < 2:
        return [1] * N

    # Load the library
    libclust = cdll.LoadLibrary("src/libclust.so")

    # Convert Python to C variables
    c_N = c_int(N)
    c_mzs = (c_double * N)(*mzs)
    d = pdist(mzs[:, np.newaxis])
    c_d = (c_double * len(d))(*d)
    c_eppm = c_double(eppm / 1e6)
    c_eabs = c_double(eabs)
    g = np.zeros(N, dtype=int)
    c_grouping = (c_int * N)(*g)

    # Calculate hierarchical peak clustering
    libclust.R_mzClust_hclust(byref(c_mzs), byref(c_N), byref(c_d), byref(c_grouping), byref(c_eppm), byref(c_eabs))
    grouping = {}
    for i in range(N):
        try:
            grouping[c_grouping[i]].append(i)
        except KeyError:
            grouping[c_grouping[i]] = [i]

    return grouping


if __name__ == "__main__":
    msfn = "/run/media/bach/EVO500GB/data/MassBank/Chubu_Univ/UT001973.txt"

    # Read all meta information from the MS-file
    spec = MBSpectrum(msfn)

    print("hey")