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
import re
import logging
import numpy as np
import os
import sqlite3
import pandas as pd

from ctypes import c_double, c_int, byref
from typing import Dict, List, Optional, Union
from scipy.spatial.distance import pdist
from zlib import crc32

import massbank2db.db
from massbank2db import HCLUST_LIB
from massbank2db.parser import get_meta_regex, get_AC_regex, get_CH_regex, get_ms_regex
from massbank2db.utils import named_row_factory
from massbank2db.exceptions import UnsupportedPrecursorType

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
        return zip(self._mz, self._int)

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

    def update_molecule_structure_information_using_pubchem(
            self, pc_dbfn: str, ids: Optional[Union[str, List[str]]] = None) -> bool:
        """
        Information of the molecular structure are extracted from PubChem. We use the CID (if provided) or the InChIKey
        provided in the Massbank file to find the compound in PubChem. The motivation is, that we have SMILES
        information coming from a single source and therefore being consistent.

        :param pc_dbfn: string, path to the local PubChem DB stored in an SQLite file

        :param ids:

        :return: boolean, indicating whether an update could be performed.
        """
        if ids is None:
            ids = ["inchi", "inchikey"]

        ids = np.atleast_1d(ids)

        # Open a connection to a local PubChemDB
        db_conn = sqlite3.connect("file:" + pc_dbfn + "?mode=ro", uri=True)  # open read-only
        db_conn.row_factory = named_row_factory

        for _id in ids:
            # Get the ID values to query information from the local PubChem DB
            id_value = self.get(_id)

            if id_value is None:
                LOGGER.info("({}) ID={} is not provided in the Massbank entry.".format(self.get("accession"), _id))
                continue

            # Fetch information from local DB (if multiple matches --> the one with lowest CID is returned)
            res = db_conn.execute(
                "SELECT cid AS pubchem_id, InChI AS inchi, InChIKey AS inchikey, SMILES_CAN AS smiles_can,"
                "       SMILES_ISO AS smiles_iso, exact_mass, molecular_formula, monoisotopic_mass, xlogp3"
                "    FROM compounds"
                "    WHERE %s is ?"
                "    ORDER BY cid ASC" % _id, (id_value, )).fetchone()

            # Update the compound with the PubChem information of found
            if res:
                for k, v in res.items():
                    self.set(k, v)

                return True

        db_conn.close()

        LOGGER.info("({}) Could not find any compound in the local PubChem DB.".format(self.get("accession")))

        return False

    def to_sirius_format(self, add_gt_molecular_formula: bool = False,
                         molecular_candidates: Optional[pd.DataFrame] = None, **kwargs) -> Dict[str, str]:
        """


        # Note: Based on [1,3,2] I would say we need to provide the 'monoisotopic_mass' to SIRIUS
        # [1] https://boecker-lab.github.io/docs.sirius.github.io/prerequisites/#monoisotopic-masses
        # [2] https://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_sdtags.pdf
        # [3] https://www.researchgate.net/post/Molecular-weight-or-exact-mass-in-LC-MS

        :param add_gt_molecular_formula:
        :param molecular_candidates:
        :param kwargs:
        :return:
        """
        output = {}
        ms_fn = self.get("accession") + ".ms"
        cand_fn = self.get("accession") + ".tsv"

        # ====================
        # Add meta-information
        # ====================
        sirius_string = ""
        sirius_string += ">compound %s\n" % self.get("accession")
        sirius_string += ">parentmass %s\n" % self.get("precursor_mz")
        sirius_string += ">ionization %s\n" % self.get("precursor_type")

        # ----------------------------------------------------------------------
        # We can provide the ground truth molecular formula to SIRIUS if desired
        # ----------------------------------------------------------------------
        _formula = "formula %s\n" % self.get("molecular_formula")
        if add_gt_molecular_formula:
            sirius_string += (">" + _formula)
        else:
            sirius_string += ("#" + _formula)  # will be ignored by SIRIUS

        if self.get("retention_time"):
            sirius_string += ">rt %f\n" % (self.get("retention_time") * 60)  # SIRIUS expects the RTs in seconds

        # -----------------------------------------
        # Try to determine the right SIRIUS profile
        # -----------------------------------------
        if "orbitrap" in self.get("instrument").lower():
            sirius_string += ">profile orbitrap\n"
        elif "qtof" in self.get("instrument_type").lower():
            sirius_string += ">profile qtof\n"
        else:
            sirius_string += ">profile default\n"

        # ----------------------------------------
        # Structure and external links as comments
        # ----------------------------------------
        sirius_string += "#smiles_iso %s\n" % self.get("smiles_iso")
        sirius_string += "#smiles_can %s\n" % self.get("smiles_can")
        sirius_string += "#inchikey %s\n" % self.get("inchikey")
        sirius_string += "#pubchem_id %s\n" % self.get("pubchem_id")
        if self.get("original_accessions") is not None:
            sirius_string += "#mb_accession %s\n" % ",".join(self.get("original_accessions"))
        else:
            sirius_string += "#mb_accession %s\n" % self.get("accession")

        # -------------------------
        # Add additional parameters
        # -------------------------
        for k, v in kwargs.items():
            sirius_string += ">{} {}\n".format(k, v)

        sirius_string += "\n"

        # =======================
        # Add fragmentation peaks
        # =======================
        if self.get("merged_peak_list", True):
            sirius_string += ">ms2merged\n"
            sirius_string += "\n".join("%f %f" % (__mz, __int) for __mz, __int in self.get_peak_list_as_tuples())
            sirius_string += "\n"
        else:
            for idx, ce in enumerate(self.get("collision_energy")):
                # sirius_string += ">collision %s\n" % ce  # TODO: Specifications like "Ramp 20.3-34eV" are not handled by SIRIUS
                sirius_string += ">ms2peaks\n"
                sirius_string += "\n".join("%f %f" % (__mz, __int)
                                           for __mz, __int in zip(self.get_mz()[idx], self.get_int()[idx]))
                sirius_string += "\n\n"

        output[ms_fn] = sirius_string

        # =================
        # Handle candidates
        # =================
        if molecular_candidates is None:
            output[cand_fn] = None
        else:
            output[cand_fn] = massbank2db.db.MassbankDB.candidates_to_sirius_format(molecular_candidates)

        return output

    def to_metfrag_format(self, molecular_candidates: Optional[pd.DataFrame] = None, **kwargs):
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
            raise UnsupportedPrecursorType("The precursor-type '%s' is not supported by MetFrag." % err)

        if not self.get("merged_peak_list", True):
            raise ValueError("MetFrag requires merged peak lists, if multiple fragmentation spectra should be processed"
                             " at ones.")

        peak_list_fn = self.get("accession") + ".peaks"
        config_fn = self.get("accession") + ".conf"

        # Peak list: tab-separated list --> mz\tint\n
        output = {
            peak_list_fn: "\n".join(["%f\t%f" % (mz, intensity) for mz, intensity in zip(self._mz, self._int)]),
        }

        # Handle candidate set
        if molecular_candidates is None:
            local_database_path = kwargs["LocalDatabasePath"]
        else:
            cands_fn = self.get("accession") + ".cands"
            local_database_path = os.path.join(kwargs["LocalDatabasePath"], cands_fn)
            output[cands_fn] = massbank2db.db.MassbankDB.candidates_to_metfrag_format(molecular_candidates)

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
                "NeutralPrecursorMass=%s" % self.get("monoisotopic_mass"),
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

                    if v[1] in ["min", "sec"]:
                        meta_info_out["retention_time_unit"] = v[1]
                    elif v[1] == "":
                        meta_info_out["retention_time_unit"] = None
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
    def merge_spectra(spectra, rt_agg_fun=np.mean, merge_peak_lists=True, normalize_peaks_before_merge=True,
                      eppm=5, eabs=0.001):
        """
        Merge the information from a list of spectra. If desired, their peak lists, e.g. belonging to different
        collision energies, are merged as well using hierarchical clustering. The meta information is merged as follows:

            - If the meta-information is equal for all spectra, then the output spectrum contains only a single info
            - If the meta-information is different for the spectra, a list of values is kept

        E.g.:

            spec_1: "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N"
            spec_2: "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N" => spec_merged "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N"
            spec_3: "inchikey" = "OUSYWCQYMPDAEO-UHFFFAOYSA-N"

            spec_1: "ce" = 10ev
            spec_1: "ce" = 20ev => spec_merged "ce" = [10ev, 20ev, 40ev]
            spec_1: "ce" = 40ev

        If the spectrum contains retention time (RT) information, than the RTs are either aggregated using the
        'rt_agg_fun' or simply concatenated if 'rt_agg_fun=None'.

        - TODO: add reference to the clustering algorithm

        :param spectra: list of MBSpectrum objects, MS/MS spectra to merge.

        :param rt_agg_fun: function, to aggregate the retention time information, if provided in the meta information of
            the spectra. The function should take in an iterable or numpy.ndarray and output a scalar

        :param merge_peak_lists: boolean, indicating whether the peak lists (i.e. fragmentation spectra) of the spectra
            should be merged using hierarchical clustering.

        :param normalize_peaks_before_merge: boolean, indicating whether the peak intensities of the individual peak
            lists to merge should be normalized (maximum value is scaled to 1.0) before merging.

        :param eppm: scalar, relative error in ppm

        :param eabs: scalar, absolute error in Da

        :return: MBSpectrum, merged spectrum
        """
        if isinstance(spectra, MBSpectrum):
            spectra = [spectra]
        elif len(spectra) == 0:
            raise ValueError("An empty spectra list cannot be merged.")

        # =====================
        # Handle the peak lists
        # =====================
        mzs_out = []
        ints_out = []

        if merge_peak_lists:
            # ----------------------------------
            # Merge them into a single peak list
            # ----------------------------------
            mzs = []
            intensities = []

            for spectrum in spectra:
                mzs += spectrum.get_mz()

                if normalize_peaks_before_merge:
                    # Work with normalized intensities
                    _norm = np.max(spectrum.get_int())
                    _ints = [_int / _norm for _int in spectrum.get_int()]
                    intensities += _ints
                else:
                    intensities += spectrum.get_int()

            mzs = np.array(mzs)
            intensities = np.array(intensities)

            # Run hierarchical clustering on the spectra peaks and return peak-grouping
            grouping = mzClust_hclust(mzs, eppm / 1e6, eabs)

            # Aggregate the peak-clusters as described in [1]: mean mass-per-charge and maximum intensities
            for peaks in grouping.values():
                mzs_out.append(np.mean(mzs[peaks]))
                ints_out.append(np.max(intensities[peaks]))
            mzs_out = np.array(mzs_out)
            ints_out = np.array(ints_out)

            # Sort the peaks (mz, int) by their mass-per-charge values
            _idc_sorted = np.argsort(mzs_out)
            mzs_out = mzs_out[_idc_sorted].tolist()
            ints_out = ints_out[_idc_sorted].tolist()
        else:
            # -------------------------------
            # Store the peak lists separately
            # -------------------------------
            for spectrum in spectra:
                mzs_out.append(spectrum.get_mz())
                ints_out.append(spectrum.get_int())

        # =================================
        # Create the merged output spectrum
        # =================================
        spec_out = MBSpectrum()
        spec_out.set("merged_peak_list", merge_peak_lists)
        spec_out.set_mz(mzs_out)
        spec_out.set_int(ints_out)

        # -------------------------------------------
        # Set meta information of the output spectrum
        # -------------------------------------------
        for info in spectra[0].get_meta_information():
            _info = []
            for spectrum in spectra:
                _info.append(spectrum.get(info))

            if len(set(_info)) == 1 \
                    and info != "retention_time" and info != "accession" and info != "collision_energy":
                spec_out.set(info, _info[0])
            else:
                spec_out.set(info, _info)

        # ----------------------------------------------------------------
        # Merge retention time information (all RTs are stored in minutes)
        # ----------------------------------------------------------------
        if rt_agg_fun is not None:
            spec_out.set("retention_time", rt_agg_fun(spec_out.get("retention_time")))

        # -------------------------------------------------------------------------------
        # Create new accession ID from the accession IDs of the individual merged spectra
        # ----------------------------------------------------------------
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

        # Ensure accession id length (10 characters)
        len_pref = len(pref)
        len_hash = 10 - len_pref
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


if __name__ == "__main__":
    msfn = "/home/bach/Documents/doctoral/data/MassBank-data_bachi55/Chubu_Univ/UT001973.txt"

    # Read all meta information from the MS-file
    spec = MBSpectrum(msfn)
