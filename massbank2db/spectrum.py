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
import os
import sqlite3

from ctypes import c_double, c_int, byref
from typing import Dict, List
from scipy.spatial.distance import pdist
from zlib import crc32

import massbank2db.db
from massbank2db import HCLUST_LIB
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

    def update_molecule_structure_information_using_pubchem(self, pc_dbfn):
        """
        Information of the molecular structure are extracted from PubChem. We use the CID (if provided) or the InChIKey
        provided in the Massbank file to find the compound in PubChem. The motivation is, that we have SMILES
        information coming from a single source and therefore being consistent.

        :param pc_dbfn: string, path to the local PubChem DB stored in an SQLite file

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
        # Open a connection to a local PubChemDB if provided
        db_conn = sqlite3.connect("file:" + pc_dbfn + "?mode=ro", uri=True)  # open read-only

        # TODO: Update also molecular weight. For that we need to rebuild the local PubChem DB.
        # TODO: Add the XLogP3 information to the spectrum
        rows = db_conn.execute("SELECT cid, InChI, InChIKey, SMILES_CAN, SMILES_ISO, exact_mass, molecular_formula "
                               "    FROM compounds"
                               "    WHERE %s is ? "
                               "    ORDER BY cid ASC" % id_type, (id, )).fetchall()

        db_conn.close()

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

    def _to_metfrag_format(self, cands=None, **kwargs):
        # All supported MetFrag precursor (ion) types: https://ipb-halle.github.io/MetFrag/projects/metfragcl/
        # Apply molecular formula simplifications, e.g. CH3COO --> C2H3O2 or CH3OH --> CH4O, source for that
        # https://docs.google.com/spreadsheets/d/1r4dPw1shIEy_W2BkfgPsihinwg-Nah654VlNTn8Gxo0/edit#gid=0
        metfrag_prec_type2mode = {
            # Positive mode
            "[M+H]+": (1, True),
            "[M+NH4]+": (18, True),
            "[M+Na]+": (23, True),
            "[M+K]+": (39, True),
            "[M+CH3OH+H]+": (33, True), "[M+CH4O+H]+": (33, True),
            "[M+ACN+H]+": (42, True), "[M+C2H3N+H]+": (42, True),
            "[M+ACN+Na]+": (64, True), "[M+C2H3N+Na]+": (64, True),
            "[M+2ACN+H]+": (83, True), "[M+C4H6N2+H]+": (83, True),
            "[M]+": (0, True),  # intrinsically charged
            # Negative mode
            "[M-H]-": (-1, False),
            "[M+Cl]-": (35, False),
            "[M+HCOO]-": (45, False), "[M+CHO2]-": (45, False),
            "[M+CH3COO]-": (59, False), "[M+C2H3O2]-": (59, False),
            "[M]-": (0, False)  # intrinsically charged
        }
        try:
            precursor_ion_mode, is_positive_mode = metfrag_prec_type2mode[self.get("precursor_type")]
        except KeyError as err:
            raise ValueError("The ionization mode '%s' is not supported by MetFrag." % err)

        peak_list_fn = self.get("accession") + ".peaks"
        config_fn = self.get("accession") + ".conf"

        # Peak list: tab-separated list --> mz\tint\n
        output = {
            peak_list_fn: "\n".join(["%f\t%f" % (mz, intensity) for mz, intensity in zip(self._mz, self._int)]),
        }

        # Handle candidate set
        if cands is None:
            local_database_path = kwargs["LocalDatabasePath"]
        else:
            cands_fn = self.get("accession") + ".cands"
            local_database_path = os.path.join(kwargs["LocalDatabasePath"], cands_fn)
            output[cands_fn] = massbank2db.db.MassbankDB._cands_to_metfrag_format(cands)

        # TODO: MetFrag configuration
        try:
            # Sanitize score information: types and weights
            score_weights = kwargs["MetFragScoreWeights"]  # type: List[int]
            score_types = kwargs["MetFragScoreTypes"]  # type: List[str]

            if not (isinstance(score_types, list) or isinstance(score_types, np.ndarray)):
                raise ValueError("Score types must be provided as list or numpy.ndarray.")

            if not (isinstance(score_weights, list) or isinstance(score_weights, np.ndarray)):
                raise ValueError("Score weights must be provided as list or numpy.ndarray.")

            if len(score_types) != len(score_weights):
                raise ValueError("Number of score types must be equal the number of score weights. (%d != %d)"
                                 % (len(score_types), len(score_weights)))

            # Sanitize pre-processing filters
            pre_processing_filters = kwargs.get("MetFragPreProcessingCandidateFilter",
                                                ["UnconnectedCompoundFilter", "IsotopeFilter"])

            if not (isinstance(pre_processing_filters, list) or isinstance(pre_processing_filters, np.ndarray)):
                raise ValueError("Pre-processing filters must be provided as list or numpy.ndarray.")

            output[config_fn] = "\n".join([
                "LocalDatabasePath=%s" % local_database_path,
                "MaximumTreeDepth=%d" % kwargs.get("MaximumTreeDepth", 2),
                "ConsiderHydrogenShifts=%s" % kwargs.get("ConsiderHydrogenShifts", True),
                "MetFragDatabaseType=%s" % kwargs.get("MetFragDatabaseType", "LocalPSV"),
                "MetFragScoreWeights=%s" % ",".join(map(str, score_weights)),
                "MetFragPreProcessingCandidateFilter=%s" % ",".join(pre_processing_filters),
                "MetFragScoreTypes=%s" % ",".join(score_types),
                "MetFragCandidateWriter=%s" % kwargs.get("MetFragCandidateWriter", "CSV"),
                "FragmentPeakMatchAbsoluteMassDeviation=%f" % kwargs.get("FragmentPeakMatchAbsoluteMassDeviation", 0.001),
                "FragmentPeakMatchRelativeMassDeviation=%f" % kwargs.get("FragmentPeakMatchRelativeMassDeviation", 5),
                "ResultsPath=%s" % kwargs["ResultsPath"],
                "NumberThreads=%d" % kwargs["NumberThreads"],
                "PrecursorIonMode=%d" % precursor_ion_mode,
                "IsPositiveIonMode=%s" % is_positive_mode,
                "NeutralPrecursorMass=%s" % self.get("exact_mass"),
                "UseSmiles=%s" % kwargs.get("UseSmiles", False),
                "SampleName=%s" % self.get("accession"),
                "PeakListPath=%s" % os.path.join(kwargs["PeakListPath"], peak_list_fn),
            ])
        except KeyError as err:
            # TODO: Should we log this event / error?
            raise ValueError("The value of '%s' is not provided and no default is defined." % err)

        return output

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

                    if v[1] in ["min", "sec", ""]:
                        meta_info_out["retention_time_unit"] = v[1]
                    elif v[1] == "s":
                        meta_info_out["retention_time_unit"] = "sec"
                    elif v[1] == "m":
                        meta_info_out["retention_time_unit"] = "min"
                    else:
                        raise ValueError("Invalid retention time-unit: '%s'" % v[1])
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
    def merge_spectra(spectra, eppm=5, eabs=0.001, rt_agg_fun=np.min):
        """
        Combining a list of spectra. Their peaks are merged into a single spectrum using hierarchical clustering. The
        meta information is merged in the following way:

            - if the meta-information is equal for all spectra, then the output spectrum contains only a single info
            - if the meta-information is different for the spectra, a list of values is kept

        E.g.:

            spec_1: "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N"
            spec_2: "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N" => spec_merged "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N"
            spec_3: "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N"

            spec_1: "ce" = 10ev
            spec_1: "ce" = 20ev => spec_merged "ce" = [10ev, 20ev, 40ev]
            spec_1: "ce" = 40ev

        If the spectrum contains retention time (RT) information (accessible via the 'rt_key' variable), than the RTs
        are either aggregated using the 'rt_agg_fun' or simply concatenated if 'rt_agg_fun=None'.

        - TODO: add reference to the clustering algorithm

        :param spectra: list of MBSpectrum objects, MS/MS spectra to merge.

        :param eppm: scalar, relative error in ppm

        :param eabs: scalar, absolute error in Da

        :param rt_agg_fun: function, to aggregate the retention time information, if provided in the meta information of
            the spectra. The function should take in an iterable or numpy.ndarray and output a scalar

        :return: MBSpectrum, merged spectrum
        """
        # Concatenate all spectra
        mzs = []
        intensities = []
        for spectrum in spectra:
            mzs += spectrum.get_mz()

            # Work with normalized intensities
            _max_int = np.max(spectrum.get_int())
            _ints = [_int / _max_int for _int in spectrum.get_int()]
            intensities += _ints

        mzs = np.array(mzs)
        intensities = np.array(intensities)

        # Run hierarchical clustering on the spectra peaks and return peak-grouping
        grouping = mzClust_hclust(mzs, eppm / 1e6, eabs)

        # Aggregate the peak-clusters as described in [1]: mean mass-per-charge and maximum intensities
        mzs_out = []
        ints_out = []
        for peaks in grouping.values():
            mzs_out.append(np.mean(mzs[peaks]))
            ints_out.append(np.max(intensities[peaks]))
        mzs_out = np.array(mzs_out)
        ints_out = np.array(ints_out)

        # Sort the peaks (mz, int) by their mass-per-charge values
        _idc_sorted = np.argsort(mzs_out)
        mzs_out = mzs_out[_idc_sorted].tolist()
        ints_out = ints_out[_idc_sorted].tolist()

        # Create the merged output spectrum
        spec_out = MBSpectrum()
        spec_out.set_mz(mzs_out)
        spec_out.set_int(ints_out)

        # Set meta information of the output spectrum
        for info in spectra[0].get_meta_information():
            _info = []
            for spectrum in spectra:
                _info.append(spectrum.get(info))

            if len(set(_info)) == 1 and info != "retention_time" and info != "accession":
                spec_out.set(info, _info[0])
            else:
                spec_out.set(info, _info)

        # Merge retention time information
        if rt_agg_fun is not None and spec_out.get("retention_time"):
            if not isinstance(spec_out.get("retention_time_unit"), str):
                raise ValueError("Merging not possible, as retention time units are not equal for all spectra: ",
                                 spec_out.get("retention_time_unit"))

            spec_out.set("retention_time", rt_agg_fun(spec_out.get("retention_time")))

        # Create new accession ID from the accession IDs of the individual merged spectra
        spec_out.set("original_accessions", spec_out.get("accession"))
        spec_out.set("accession", MBSpectrum._get_new_accession_id(spec_out.get("accession")))

        return spec_out

    @staticmethod
    def _get_new_accession_id(accs):
        # Extract the accession prefix, e.g. AU203981 -> AU or FEL02392 -> FEL
        regex = re.compile("[A-Z]+")
        pref = [regex.match(acc)[0] for acc in accs]
        assert all([_pref == pref[0] for _pref in pref]), "All accession prefixes are assumed to be equal."
        pref = pref[0]

        # Get CRC32 hash value based on all accession ids
        hash_int = crc32("".join(accs).encode('utf-8'))

        # Ensure accession id length (8 characters)
        len_pref = len(pref)
        len_hash = 8 - len_pref
        hash_str = ("%%0%dd" % len_hash) % (hash_int % (10 ** len_hash - 1))

        # Combine the accession predix with the hash sting
        return pref + hash_str


def mzClust_hclust(mzs, eppm, eabs):
    N = len(mzs)

    if N < 2:
        # Empty spectrum or single peak: no grouping needed
        c_grouping = np.ones(N)
    else:
        # Convert Python to C variables
        c_N = c_int(N)
        c_mzs = (c_double * N)(*mzs)
        d = pdist(mzs[:, np.newaxis])
        c_d = (c_double * len(d))(*d)
        c_eppm = c_double(eppm)
        c_eabs = c_double(eabs)
        g = np.zeros(N, dtype=int)
        c_grouping = (c_int * N)(*g)

        # Calculate hierarchical peak clustering
        HCLUST_LIB.R_mzClust_hclust(byref(c_mzs), byref(c_N), byref(c_d), byref(c_grouping), byref(c_eppm),
                                    byref(c_eabs))

    grouping = {}
    for i in range(N):
        try:
            grouping[c_grouping[i]].append(i)
        except KeyError:
            grouping[c_grouping[i]] = [i]

    return grouping

# OLD CODE
# Apply filter to the each spectrum as described in [1]
# 1) Keep only spectra, that contain at least one fragment peak _other_ than the precursor
# _mzs = spectrum.get_mz()
# if len(_mzs) < 2:
#     print("skip:", spectrum.get("accession"))
#     continue
# if all([_within_error_window_around_peak(float(spectrum.get("precursor_mz")), _mz, eppm, eabs)
#         for _mz in _mzs]):
#     print("skip:", spectrum.get("accession"))
#     continue

# def _within_error_window_around_peak(mz_center, mz, eppm, eabs):
#     lb = mz_center - mz_center * eppm - eabs
#     up = mz_center + mz_center * eppm + eabs
#
#     if (lb <= mz) and (mz <= up):
#         return True
#     else:
#         return False


if __name__ == "__main__":
    msfn = "/run/media/bach/EVO500GB/data/MassBank/Chubu_Univ/UT001973.txt"

    # Read all meta information from the MS-file
    spec = MBSpectrum(msfn)
