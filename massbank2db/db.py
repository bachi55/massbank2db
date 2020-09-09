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

import sqlite3
import os
import glob
import re
import logging
import pandas as pd
import numpy as np

from hashlib import sha1

from massbank2db.spectrum import MBSpectrum

# Setup the Logger
LOGGER = logging.getLogger(__name__)
CH = logging.StreamHandler()
FORMATTER = logging.Formatter('[%(levelname)s] %(name)s : %(message)s')
CH.setFormatter(FORMATTER)
LOGGER.addHandler(CH)


def create_db(file_pth):
    """
    Create an SQLite DB to store the Massbank data.

    :param file_pth: string, DB file path.

    Example:
        >>> from massbank2db.db import create_db
        >>> db_pth = 'library.db'
        >>> create_db(file_pth=db_pth)
    """
    with sqlite3.connect(file_pth) as conn:
        for table_name in ["spectra_candidates", "spectra_peaks", "spectra_raw_rts", "spectra_meta", "datasets",
                           "molecules"]:
            conn.execute("DROP TABLE IF EXISTS %s" % table_name)

        # Molecules Table
        conn.execute(
            "CREATE TABLE molecules( \
                 cid               INTEGER PRIMARY KEY NOT NULL, \
                 inchi             VARCHAR NOT NULL, \
                 inchikey          VARCHAR NOT NULL, \
                 inchikey1         VARCHAR NOT NULL, \
                 inchikey2         VARCHAR NOT NULL, \
                 smiles_iso        VARCHAR NOT NULL, \
                 smiles_can        VARCHAR NOT NULL, \
                 exact_mass        FLOAT NOT NULL, \
                 molecular_formula VARCHAR NOT NULL)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey_index ON molecules(inchikey)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey1_index ON molecules(inchikey1)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey2_index ON molecules(inchikey2)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_exact_mass_index ON molecules(exact_mass)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_mf_index ON molecules(molecular_formula)")

        # Datasets Meta-information table
        #   primary key: Accession prefix + some running id
        #                E.g. AU_001, AU_002
        conn.execute(
            "CREATE TABLE datasets( \
                name                VARCHAR PRIMARY KEY NOT NULL, \
                contributor         VARCHAR NOT NULL, \
                ion_mode            VARCHAR NOT NULL, \
                num_spectra         INTEGER, \
                num_compounds       INTEGER, \
                copyright           VARCHAR, \
                license             VARCHAR, \
                column_name         VARCHAR NOT NULL, \
                column_type         VARCHAR NOT NULL, \
                column_temperature  FLOAT, \
                flow_gradient       VARCHAR NOT NULL, \
                flow_rate           FLOAT, \
                solvent_A           VARCHAR NOT NULL, \
                solvent_B           VARCHAR NOT NULL, \
                solvent             VARCHAR, \
                instrument_type     VARCHAR NOT NULL, \
                instrument          VARCHAR NOT NULL, \
                column_dead_time    FLOAT)"
        )

        # Spectra Meta-information table
        # Note:
        #   "[...] if the foreign key column in the track table is NULL, then no corresponding entry in the artist
        #   table is required."
        #   Source: https://sqlite.org/foreignkeys.html#fk_basics
        conn.execute(
            "CREATE TABLE spectra_meta( \
                accession           VARCHAR PRIMARY KEY NOT NULL, \
                dataset             VARCHAR NOT NULL, \
                record_title        VARCHAR NOT NULL, \
                molecule            INTEGER NOT NULL, \
                precursor_mz        FLOAT NOT NULL, \
                precursor_type      FLOAT NOT NULL, \
                collision_energy    FLOAT, \
                ms_type             VARCHAR NOT NULL, \
                resolution          FLOAT, \
                fragmentation_mode  VARCHAR, \
                spectra_group       VARCHAR, \
             FOREIGN KEY(molecule)      REFERENCES molecules(cid),\
             FOREIGN KEY(dataset)       REFERENCES datasets(name) ON DELETE CASCADE)"
        )
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_meta_dataset_index ON spectra_meta(dataset)")
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_meta_spectra_group_index ON spectra_meta(spectra_group)")

        # Spectra Retention Time information table
        conn.execute(
            "CREATE TABLE spectra_raw_rts( \
                spectrum            VARCHAR NOT NULL, \
                retention_time      FLOAT NOT NULL, \
                retention_time_unit VARCHAR DEFAULT 'min', \
             FOREIGN KEY(spectrum) REFERENCES spectra_meta(accession) ON DELETE CASCADE)"
        )
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_raw_rts_spectrum_index ON spectra_raw_rts(spectrum)")
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_raw_rts_retention_time_index ON spectra_raw_rts(retention_time)")

        # Spectra Peak Information table
        conn.execute(
            "CREATE TABLE spectra_peaks( \
                spectrum  VARCHAR NOT NULL, \
                mz        FLOAT NOT NULL, \
                intensity FLOAT NOT NULL, \
             FOREIGN KEY(spectrum) REFERENCES spectra_meta(accession) ON DELETE CASCADE)"
        )
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_peaks_spectrum_index ON spectra_peaks(spectrum)")

        # Spectra candidate information
        conn.execute(
            "CREATE TABLE spectra_candidates( \
                spectrum      VARCHAR NOT NULL, \
                candidate     INTEGER NOT NULL, \
                mf_gt_equal   INTEGER NOT NULL, \
                ppm_diff_gt   FLOAT NOT NULL, \
             FOREIGN KEY(spectrum)  REFERENCES spectra_meta(accession)  ON DELETE CASCADE, \
             FOREIGN KEY(candidate) REFERENCES molecules(cid)           ON DELETE CASCADE, \
             PRIMARY KEY(spectrum, candidate))")
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_candidates_spectrum_index ON spectra_candidates(spectrum)")
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_candidates_candidate_index ON spectra_candidates(candidate)")
        conn.execute("CREATE INDEX IF NOT EXISTS spectra_candidates_ppm_diff_index ON spectra_candidates(ppm_diff_gt)")


def get_temporal_database(file_pth=":memory:") -> sqlite3.Connection:
    conn = sqlite3.connect(file_pth)
    conn.execute("DROP TABLE IF EXISTS information")
    conn.execute("CREATE TABLE information ("
                 "  accession           VARCHAR NOT NULL PRIMARY KEY,"
                 "  contributor         VARCHAR NOT NULL,"
                 "  accession_prefix    VARCHAR NOT NULL,"
                 "  instrument_type     VARCHAR NOT NULL,"
                 "  instrument          VARCHAR NOT NULL,"
                 "  ion_mode            VARCHAR NOT NULL,"
                 "  column_name         VARCHAR,"
                 "  flow_gradient       VARCHAR,"
                 "  flow_rate           VARCHAR,"
                 "  solvent_A           VARCHAR,"
                 "  solvent_B           VARCHAR,"
                 "  solvent             VARCHAR,"
                 "  column_temperature  VARCHAR,"
                 "  inchikey            VARCHAR NOT NULL,"
                 "  ms_level            VARCHAR,"       
                 "  retention_time      FLOAT)")

    return conn


class MassbankDB(object):
    def __init__(self, mb_dbfn, only_with_rt=True, only_ms2=True, use_pubchem_structure_info=True,
                 exclude_deprecated=True, min_number_of_unique_compounds_per_dataset=50, pc_dbfn=None):
        """
        
        :param mb_dbfn:
        :param only_with_rt:
        :param only_ms2:
        :param use_pubchem_structure_info:
        :param min_number_of_unique_compounds_per_dataset:
        :param pc_dbfn:
        """
        self.__mb_conn = sqlite3.connect(mb_dbfn)

        # Open a connection to a local PubChemDB if provided
        if pc_dbfn is not None:
            self.__pc_conn = sqlite3.connect("file:" + pc_dbfn + "?mode=ro", uri=True)  # open read-only
        else:
            self.__pc_conn = None

        # Different filters to exclude entries from the database
        self.__only_with_rt = only_with_rt
        self.__only_ms2 = only_ms2
        self.__use_pubchem_structure_info = use_pubchem_structure_info
        self.__min_number_of_unique_compounds_per_dataset = min_number_of_unique_compounds_per_dataset
        self.__exclude_deprecated = exclude_deprecated

    def __enter__(self):
        """
        When entering the context manager
        """
        # Enable Foreign Key Support (see also: https://sqlite.org/foreignkeys.html#fk_enable)
        self.__mb_conn.execute("PRAGMA foreign_keys = ON")

        return self

    def __exit__(self, ext_type, exc_value, traceback):
        """
        When leaving the context manager
        """
        # Rollback Massbank database in case of an error
        if isinstance(exc_value, Exception):
            self.__mb_conn.rollback()

        # Close connection to the Massbank database
        self.__mb_conn.close()

        # Close connection to the PubChem database
        if self.__pc_conn:
            self.__pc_conn.close()

    def insert_dataset(self, accession_prefix, contributor, base_path):
        """
        Insert all accession of the specified contributor with specified prefix to the database. Each contributor and
        each corresponding accession prefix represents a separate dataset. Furthermore, different chromatographic (LC)
        configurations, MS ionization modes and MS instruments are split into sub-datasets.

        Example:
            (AU, Athens_Univ)
                -> AU_0001  ~ positive, LC configuration 1, MS instrument 1
                -> AU_0002  ~ positive, LC configuration 2, MS instrument 1
                -> AU_1001  ~ negative, LC configuration 3, MS instrument 1
                -> AU_1002  ~ negative, LC configuration 4, MS instrument 2
                -> ...

        The LC configuration, MS ionization mode and MS instrument is stored as meta data in the dataset table.

        :param accession_prefix:
        :param contributor:
        :param base_path:
        :return:
        """
        def _add_spec_to_filter_db(spec):
            # Specify the solvent information
            solvent_A = spec.get("solvent_A")
            solvent_B = spec.get("solvent_B")
            if solvent_A and solvent_B:
                solvent = None
            else:
                solvent = spec.get("solvent")

            new_row = (acc, contributor, accession_prefix, spec.get("instrument_type"), spec.get("instrument"),
                       spec.get("ion_mode"), spec.get("column_name"), spec.get("flow_gradient"), spec.get("flow_rate"),
                       solvent_A, solvent_B, solvent, spec.get("column_temperature"), spec.get("inchikey"),
                       spec.get("ms_type"), spec.get("retention_time"))

            filter_db_conn.execute("INSERT INTO information VALUES(%s)" % self._get_db_value_placeholders(len(new_row)),
                                   new_row)

        # ================================================================
        # Build up a temporary database to filter and separate accessions:
        # ================================================================
        # - group accessions by their MS and LC configuration --> each combination becomes a separate dataset
        # - remove groups of accessions, where we have less than a certain number
        #   (see 'min_number_of_unique_compounds_per_dataset')
        filter_db_conn = get_temporal_database()  # in memory
        specs = {}
        with filter_db_conn:
            for msfn in glob.iglob(os.path.join(base_path, contributor, accession_prefix + "[0-9]*.txt")):
                acc = os.path.basename(msfn).split(".")[0]
                specs[acc] = MBSpectrum(msfn)
                assert specs[acc].get("accession") == acc

                if specs[acc].get("inchikey") is None:
                    LOGGER.info("(%s) No compound structure information, i.e. no InChIKey." % acc)
                    continue

                if specs[acc].get("deprecated") is not None:
                    LOGGER.info("(%s) Deprecated entry: '%s'" % (acc, specs[acc].get("deprecated")))
                    continue

                _add_spec_to_filter_db(specs[acc])

        # =====================
        # Group the accessions:
        # =====================
        _stmt = ["SELECT COUNT(DISTINCT inchikey), GROUP_CONCAT(accession) FROM information"]
        if self.__only_with_rt and self.__only_ms2:
            _stmt.append("WHERE retention_time IS NOT NULL AND ms_level IS 'MS2'")
        elif self.__only_with_rt:
            _stmt.append("WHERE retention_time IS NOT NULL")
        elif self.__only_ms2:
            _stmt.append("WHERE ms_level IS 'MS2'")
        _stmt.append("GROUP BY ion_mode, instrument, instrument_type, column_name, column_temperature, flow_gradient, "
                     "  flow_rate, solvent_A, solvent_B")

        cur = filter_db_conn.execute(" ".join(_stmt))

        # ============================================
        # Include the accessions as separate datasets:
        # ============================================
        acc_pref_idx = 0
        for row in cur:
            if row[0] < self.__min_number_of_unique_compounds_per_dataset:
                continue

            dataset_identifier = "%s_%03d" % (accession_prefix, acc_pref_idx)

            for acc in row[1].split(","):
                # -----------------------------------------------------------------------------
                # Update the molecule structure information in the spectrum file using PubChem:
                # -----------------------------------------------------------------------------
                if self.__use_pubchem_structure_info and \
                        not specs[acc].update_molecule_structure_information_using_pubchem(self.__pc_conn):
                    continue

                # -----------------------------------
                # Insert the spectrum to the database
                # -----------------------------------
                try:
                    with self.__mb_conn:
                        self.insert_spectrum(dataset_identifier, contributor, specs[acc])
                except sqlite3.IntegrityError as err:
                    LOGGER.error("({}) SQLite integrity error: {}".format(acc, err))
                    continue

            # ----------------------------------------------------
            # Determine the number of spectra and unique compounds
            # ----------------------------------------------------
            with self.__mb_conn:
                _num_spec, _num_cmp = self.__mb_conn.execute("SELECT COUNT(accession), COUNT(distinct molecule) "
                                                             "   FROM spectra_meta "
                                                             "   WHERE dataset IS ?", (dataset_identifier,)).fetchall()[0]
                self.__mb_conn.execute("UPDATE datasets "
                                       "   SET num_spectra = ?, num_compounds = ?"
                                       "   WHERE name IS ?", (_num_spec, _num_cmp, dataset_identifier))

            # -----------------------------
            # Group the spectra for merging
            # -----------------------------
            with self.__mb_conn:
                self.group_spectra()

            acc_pref_idx += 1

    def insert_spectrum(self, dataset_identifier: str, contributor: str, spectrum: MBSpectrum):
        """

        :param spectrum:
        :return:
        """
        # ===========================
        # Insert Molecule Information
        # ===========================
        self.__mb_conn.execute("INSERT OR IGNORE INTO molecules VALUES (%s)" % self._get_db_value_placeholders(9),
                               (
                                   spectrum.get("pubchem_id"),
                                   spectrum.get("inchi"),
                                   spectrum.get("inchikey"),
                                   spectrum.get("inchikey").split("-")[0],
                                   spectrum.get("inchikey").split("-")[1],
                                   spectrum.get("smiles_iso"),
                                   spectrum.get("smiles_can"),
                                   spectrum.get("exact_mass"),
                                   spectrum.get("molecular_formula")
                               ))

        # ==========================
        # Insert Dataset Information
        # ==========================
        ion_mode = spectrum.get("ion_mode")
        if not ion_mode:
            # have to do special check for ionization mode (as sometimes gets missed)
            m = re.search("^\[.*\](\-|\+)", spectrum.get("precursor_type"), re.IGNORECASE)
            if m:
                polarity = m.group(1).strip()
                if polarity == "+":
                    ion_mode = "positive"
                elif polarity == "-":
                    ion_mode = "negative"
                else:
                    raise RuntimeError("Ups")
        else:
            ion_mode = ion_mode.lower()
        self.__mb_conn.execute("INSERT OR IGNORE INTO datasets VALUES (%s)" % self._get_db_value_placeholders(18),
                               (
                                   dataset_identifier,
                                   contributor,
                                   ion_mode,
                                   0,
                                   0,
                                   spectrum.get("copyright"),
                                   spectrum.get("license"),
                                   spectrum.get("column_name"),
                                   "",  # FIXME: column type, e.g. RP or HILIC, is not specified in the Massbank file
                                   spectrum.get("column_temperature"),
                                   spectrum.get("flow_gradient"),
                                   spectrum.get("flow_rate"),
                                   spectrum.get("solvent_A"),
                                   spectrum.get("solvent_B"),
                                   spectrum.get("solvent"),
                                   spectrum.get("instrument_type"),
                                   spectrum.get("instrument"),
                                   None  # TODO: the column dead-time needs to be estimated from the data
                               ))

        # ===============================
        # Insert Spectra Meta Information
        # ===============================
        self.__mb_conn.execute("INSERT INTO spectra_meta VALUES (%s)" % self._get_db_value_placeholders(11),
                               (
                                   spectrum.get("accession"),
                                   dataset_identifier,
                                   spectrum.get("record_title"),
                                   spectrum.get("pubchem_id"),
                                   spectrum.get("precursor_mz"),
                                   spectrum.get("precursor_type"),
                                   spectrum.get("collision_energy"),
                                   spectrum.get("ms_type"),
                                   spectrum.get("resolution"),
                                   spectrum.get("fragmentation_mode"),
                                   None
                               ))

        # ====================
        # Insert Spectra Peaks
        # ====================
        _mz, _int = spectrum.get_mz(), spectrum.get_int()
        _n_peaks = len(_mz)
        _acc = [spectrum.get("accession")] * _n_peaks
        self.__mb_conn.executemany("INSERT INTO spectra_peaks VALUES(%s)" % self._get_db_value_placeholders(3),
                                   zip(_acc, _mz, _int))

        # =====================
        # Insert Retention Time
        # =====================
        if spectrum.get("retention_time_unit") == '':
            if spectrum.get("retention_time") > 100:
                LOGGER.warning("(%s) Retention time unit default (=min) does not look reasonable: rt=%.3f"
                               % (spectrum.get("accession"), spectrum.get("retention_time")))

            self.__mb_conn.execute("INSERT INTO spectra_raw_rts (spectrum, retention_time) VALUES(?, ?)",
                                   (
                                    spectrum.get("accession"),
                                    spectrum.get("retention_time"),
                                   ))
        else:
            self.__mb_conn.execute("INSERT INTO spectra_raw_rts VALUES(?, ?, ?)",
                                   (
                                    spectrum.get("accession"),
                                    spectrum.get("retention_time"),
                                    spectrum.get("retention_time_unit")
                                   ))

    def group_spectra(self):
        """
        Group spectra such that each group represents a set of spectra to merged prior analysis, e.g. using MetFrag or
        CSI:FingerID.
        """
        rows = self.__mb_conn.execute("SELECT dataset, \"('\" || GROUP_CONCAT(accession, \"','\") || \"')\" "
                                      "   FROM spectra_meta GROUP BY dataset, molecule, precursor_mz, precursor_type,"
                                      "   resolution, fragmentation_mode")

        for ds, l_accs in rows:
            acc_grp_id = ds + "_" + sha1(l_accs.encode('utf-8')).hexdigest()[:8]   # e.g. AU_001_3a1fd8de
            self.__mb_conn.execute("UPDATE spectra_meta SET spectra_group = ? "
                                   "    WHERE accession IN %s" % l_accs, (acc_grp_id, ))

    def iter_spectra(self, dataset, grouped=True, return_candidates=False, ppm=5):
        if return_candidates and not self.__pc_conn:
            raise ValueError("No local PubChem DB. Candidates cannot be queried.")

        if grouped:
            # rows = self.__mb_conn.execute("SELECT GROUP_CONCAT(accession), dataset, molecule, precursor_mz,"
            #                               "       precursor_type, GROUP_CONCAT(collision_energy),"
            #                               "       GROUP_CONCAT(ms_type), resolution, fragmentation_mode "
            #                               "   FROM spectra_meta "
            #                               "   WHERE dataset IS ? "
            #                               "   GROUP BY dataset, molecule, precursor_mz, precursor_type, resolution,"
            #                               "            fragmentation_mode", (dataset, ))

            rows = self.__mb_conn.execute("SELECT GROUP_CONCAT(accession), dataset, molecule, precursor_mz,"
                                          "       precursor_type, GROUP_CONCAT(collision_energy),"
                                          "       GROUP_CONCAT(ms_type), GROUP_CONCAT(resolution), fragmentation_mode "
                                          "   FROM spectra_meta "
                                          "   WHERE dataset IS ? "
                                          "   GROUP BY dataset, molecule, precursor_mz, precursor_type,"
                                          "            fragmentation_mode", (dataset,))
        else:
            rows = self.__mb_conn.execute("SELECT accession, dataset, molecule, precursor_mz, precursor_type, "
                                          "       collision_energy, ms_type, resolution, fragmentation_mode "
                                          "   FROM spectra_meta "
                                          "   WHERE dataset IS ?", (dataset,))

        for row in rows:
            specs = []

            # -------------------------
            # Load compound information
            # -------------------------
            # Note: This is equal for all spectra, if grouped, as the compound is a grouping criteria.
            mol = self.__mb_conn.execute("SELECT * FROM molecules WHERE cid = ?", (row[2],)).fetchall()[0]

            # --------------------------
            # Load candidate information
            # --------------------------
            if return_candidates is None:
                cands = None
            elif return_candidates == "mf":
                cands = self.__pc_conn.execute("SELECT cid, smiles_iso FROM compounds WHERE molecular_formula IS ?",
                                               (mol[8],)).fetchall()
            elif return_candidates == "mz":
                min_exact_mass, max_exact_mass = self._get_ppm_window(mol[7], ppm)
                cands = self.__pc_conn.execute("SELECT cid, smiles_iso FROM compounds WHERE exact_mass BETWEEN ? AND ?",
                                               (min_exact_mass, max_exact_mass)).fetchall()
            else:
                raise ValueError("Invalid")

            # -----------------------------------------------------
            # Create a Spectrum object for all spectra in the group
            # -----------------------------------------------------
            for acc, ce, ms_type in zip(row[0].split(","), row[5].split(","), row[6].split(",")):
                specs.append(MBSpectrum())

                # --------------------------
                # Data about the acquisition
                # --------------------------
                specs[-1].set("accession", acc)
                specs[-1].set("dataset", row[1])
                specs[-1].set("precursor_mz", row[3])
                specs[-1].set("precursor_type", row[4])
                specs[-1].set("collision_energy", ce)
                specs[-1].set("ms_type", ms_type)
                specs[-1].set("resolution", row[7])
                specs[-1].set("fragmentation_mode", row[8])

                # -------------------
                # Retention time data
                # -------------------
                rt, rt_unit = self.__mb_conn.execute("SELECT retention_time, retention_time_unit FROM spectra_raw_rts"
                                                     "   WHERE spectrum IS ?", (acc, )).fetchall()[0]
                specs[-1].set("retention_time", rt)
                specs[-1].set("retention_time_unit", rt_unit)

                # --------------------
                # Compound information
                # --------------------
                specs[-1].set("pubchem_id", mol[0])
                specs[-1].set("inchi", mol[1])
                specs[-1].set("inchikey", mol[2])
                specs[-1].set("smiles_iso", mol[5])
                specs[-1].set("smiles_can", mol[6])
                specs[-1].set("exact_mass", mol[7])
                specs[-1].set("molecular_formula", mol[8])

                # --------------
                # Spectrum peaks
                # --------------
                # TODO: Remove loop here
                peaks = self.__mb_conn.execute("SELECT mz, intensity FROM spectra_peaks "
                                               "    WHERE spectrum IS ? ORDER BY mz", (acc, ))
                specs[-1]._mz, specs[-1]._int = [], []
                for peak in peaks:
                    specs[-1]._mz.append(peak[0])
                    specs[-1]._int.append(peak[1])

            yield mol, specs, cands

    @staticmethod
    def _get_ppm_window(exact_mass, ppm):
        abs_deviation = exact_mass / 1e6 * ppm
        return exact_mass - abs_deviation, exact_mass + abs_deviation

    @staticmethod
    def _get_db_value_placeholders(n, placeholder='?'):
        """

        :param n:
        :return:
        """
        return ",".join([placeholder] * n)

    @staticmethod
    def _in_sql(li):
        """
        Concatenates a list of strings to a SQLite ready string that can be used in combination
        with the 'in' statement.

        E.g.:
            ["house", "boat", "donkey"] --> "('house', 'boat', 'donkey')"


        :param li: list of strings

        :return: SQLite ready string for 'in' statement
        """
        return "(" + ",".join(["'%s'" % li for li in np.atleast_1d(li)]) + ")"

            
if __name__ == "__main__":
    LOGGER.setLevel(logging.INFO)

    massbank_dir = "/run/media/bach/EVO500GB/data/MassBank"

    dbfn = "tests/test_DB.sqlite"
    pubchem_dbfn = "/run/media/bach/EVO500GB/data/pubchem_24-06-2019/db/pubchem.sqlite"

    create_db(dbfn)

    # Load list of datasets provided by Massbank
    mbds = pd.read_csv(os.path.join(massbank_dir, "List_of_Contributors_Prefixes_and_Projects.md"),
                       sep="|", skiprows=2, header=None) \
        .iloc[:, [1, 4]] \
        .applymap(str.strip) \
        .rename({1: "Contributor", 4: "AccPref"}, axis=1)  # type: pd.DataFrame

    for idx, row in mbds.iterrows():
        print("(%02d/%02d) %s: " % (idx + 1, len(mbds), row["Contributor"]))
        for pref in map(str.strip, row["AccPref"].split(",")):
            print(pref)
            with MassbankDB(dbfn, pc_dbfn=pubchem_dbfn) as mbdb:
                mbdb.insert_dataset(pref, row["Contributor"], massbank_dir)