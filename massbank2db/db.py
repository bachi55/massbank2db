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

import sqlite3
import os
import glob
import re
import logging
import pandas as pd
import numpy as np

from typing import Optional, List

import massbank2db.spectrum
from massbank2db.utils import get_mass_error_in_ppm

# Setup the Loggers
LOGGER = logging.getLogger(__name__)
__SH = logging.StreamHandler()
__SH.setFormatter(logging.Formatter('[%(levelname)s] %(name)s : %(message)s'))
LOGGER.addHandler(__SH)


class MassbankDB(object):
    def __init__(self, mb_dbfn, read_only=False):
        """
        :param mb_dbfn:
        """
        if read_only:
            self._mb_conn = sqlite3.connect("file:" + mb_dbfn + "?mode=ro", uri=True)
        else:
            self._mb_conn = sqlite3.connect(mb_dbfn)

    def __enter__(self):
        """
        When entering the context manager
        """
        # Enable Foreign Key Support (see also: https://sqlite.org/foreignkeys.html#fk_enable)
        self._mb_conn.execute("PRAGMA foreign_keys = ON")

        return self

    def __exit__(self, ext_type, exc_value, traceback):
        """
        When leaving the context manager
        """
        # Rollback Massbank database in case of an error
        if isinstance(exc_value, Exception):
            self._mb_conn.rollback()

        # Close connection to the Massbank database
        self.close()

    def close(self):
        self._mb_conn.close()

    def get_datasets_table(self):
        """
        Return the datasets table as pandas.DataFrame
        """
        return pd.read_sql_query("SELECT * FROM datasets", con=self._mb_conn)

    def initialize_tables(self, reset=False):
        """
        Initialize the database tables.

        :param reset: boolean, indicating whether the existing tables should be dropped.
        """
        with self._mb_conn:
            if reset:
                for table_name in ["spectra_peaks", "spectra_raw_rts", "spectra_meta", "datasets", "molecules"]:
                    self._mb_conn.execute("DROP TABLE IF EXISTS %s" % table_name)

            # Molecules Table
            self._mb_conn.execute(
                "CREATE TABLE IF NOT EXISTS molecules( \
                     cid               INTEGER PRIMARY KEY NOT NULL, \
                     inchi             VARCHAR NOT NULL, \
                     inchikey          VARCHAR NOT NULL, \
                     inchikey1         VARCHAR NOT NULL, \
                     inchikey2         VARCHAR NOT NULL, \
                     smiles_iso        VARCHAR NOT NULL, \
                     smiles_can        VARCHAR, \
                     exact_mass        FLOAT, \
                     monoisotopic_mass FLOAT NOT NULL, \
                     molecular_formula VARCHAR NOT NULL, \
                     xlogp3            FLOAT)")

            # Datasets Meta-information table
            #   primary key: Accession prefix + some running id
            #                E.g. AU_001, AU_002
            self._mb_conn.execute(
                "CREATE TABLE IF NOT EXISTS datasets( \
                    name                VARCHAR PRIMARY KEY NOT NULL, \
                    contributor         VARCHAR NOT NULL, \
                    ion_mode            VARCHAR NOT NULL, \
                    num_spectra         INTEGER, \
                    num_compounds       INTEGER, \
                    copyright           VARCHAR, \
                    license             VARCHAR, \
                    column_name         VARCHAR NOT NULL, \
                    column_type         VARCHAR, \
                    column_temperature  FLOAT, \
                    flow_gradient       VARCHAR NOT NULL, \
                    flow_rate           FLOAT, \
                    solvent_A           VARCHAR NOT NULL, \
                    solvent_B           VARCHAR NOT NULL, \
                    solvent             VARCHAR, \
                    instrument_type     VARCHAR NOT NULL, \
                    instrument          VARCHAR NOT NULL, \
                    fragmentation_mode  VARCHAR, \
                    column_dead_time    FLOAT)"
            )

            # Spectra Meta-information table
            # Note:
            #   "[...] if the foreign key column in the track table is NULL, then no corresponding entry in the artist
            #   table is required."
            #   Source: https://sqlite.org/foreignkeys.html#fk_basics
            self._mb_conn.execute(
                "CREATE TABLE IF NOT EXISTS spectra_meta( \
                    accession            VARCHAR PRIMARY KEY NOT NULL, \
                    dataset              VARCHAR NOT NULL, \
                    record_title         VARCHAR NOT NULL, \
                    molecule             INTEGER NOT NULL, \
                    precursor_mz         FLOAT NOT NULL, \
                    precursor_type       FLOAT NOT NULL, \
                    collision_energy     FLOAT, \
                    ms_type              VARCHAR NOT NULL, \
                    resolution           FLOAT, \
                    fragmentation_mode   VARCHAR, \
                    exact_mass_error_ppm FLOAT, \
                 FOREIGN KEY(molecule)      REFERENCES molecules(cid),\
                 FOREIGN KEY(dataset)       REFERENCES datasets(name) ON DELETE CASCADE)"
            )

            # Spectra Retention Time information table
            self._mb_conn.execute(
                "CREATE TABLE IF NOT EXISTS spectra_raw_rts( \
                    spectrum            VARCHAR NOT NULL, \
                    retention_time      FLOAT NOT NULL, \
                    retention_time_unit VARCHAR DEFAULT 'min', \
                 FOREIGN KEY(spectrum) REFERENCES spectra_meta(accession) ON DELETE CASCADE)"
            )

            # Spectra Peak Information table
            self._mb_conn.execute(
                "CREATE TABLE IF NOT EXISTS spectra_peaks( \
                    spectrum  VARCHAR NOT NULL, \
                    mz        FLOAT NOT NULL, \
                    intensity FLOAT NOT NULL, \
                 FOREIGN KEY(spectrum) REFERENCES spectra_meta(accession) ON DELETE CASCADE)"
            )

            # ------------------
            # Create indices
            # ------------------
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchi_index             ON molecules(inchi)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey_index          ON molecules(inchikey)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey1_index         ON molecules(inchikey1)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey2_index         ON molecules(inchikey2)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS molecules_exact_mass_index        ON molecules(exact_mass)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS molecules_monoisotopic_mass_index ON molecules(monoisotopic_mass)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS molecules_mf_index                ON molecules(molecular_formula)")

            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS spectra_meta_dataset_index           ON spectra_meta(dataset)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS spectra_raw_rts_spectrum_index       ON spectra_raw_rts(spectrum)")
            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS spectra_raw_rts_retention_time_index ON spectra_raw_rts(retention_time)")

            self._mb_conn.execute("CREATE INDEX IF NOT EXISTS spectra_peaks_spectrum_index ON spectra_peaks(spectrum)")

            # ------------------
            # Create triggers
            # ------------------
            # Remove molecules that are not referenced anymore by any spectra, e.g. after a filtering operation
            self._mb_conn.execute("CREATE TRIGGER remove_molecules_without_spectra_reference "
                                  "    AFTER DELETE "
                                  "    ON spectra_meta"
                                  "    BEGIN"
                                  "        DELETE FROM molecules WHERE cid IN ("
                                  "            SELECT cid FROM molecules "
                                  "                LEFT OUTER JOIN spectra_meta sm ON molecules.cid = sm.molecule"
                                  "                WHERE accession IS NULL"
                                  "            );"
                                  "    END")

    def insert_dataset(self, accession_prefix, contributor, base_path, only_with_rt=True, only_ms2=True,
                       use_pubchem_structure_info=True, pc_dbfn=None, exclude_deprecated=True, only_lc=True):
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
                       spec.get("ms_type"), spec.get("retention_time"), spec.get("fragmentation_mode"))

            filter_db_conn.execute("INSERT INTO information VALUES(%s)" % self._get_db_value_placeholders(len(new_row)),
                                   new_row)

        # ================================================================
        # Build up a temporary database to filter and separate accessions:
        # ================================================================
        # - group accessions by their MS and LC configuration --> each combination becomes a separate dataset
        # - remove groups of accessions, where we have less than a certain number
        #   (see 'min_number_of_unique_compounds_per_dataset')
        filter_db_conn = self._get_temporal_database()  # in memory
        specs = {}

        msfns = glob.glob(os.path.join(base_path, contributor, accession_prefix + "[0-9]*.txt"))
        LOGGER.info("Number of files: %d" % len(msfns))
        for msfn in msfns:
            acc = os.path.basename(msfn).split(".")[0]
            specs[acc] = massbank2db.spectrum.MBSpectrum(msfn)
            assert specs[acc].get("accession") == acc

            if specs[acc].get("inchikey") is None:
                LOGGER.info("(%s) No compound structure information, i.e. no InChIKey." % acc)
                continue

            if exclude_deprecated and specs[acc].get("deprecated") is not None:
                LOGGER.info("(%s) Deprecated entry: '%s'" % (acc, specs[acc].get("deprecated")))
                continue

            _add_spec_to_filter_db(specs[acc])

        # =====================
        # Group the accessions:
        # =====================
        _stmt = ["SELECT COUNT(DISTINCT inchikey), GROUP_CONCAT(accession) FROM information"]

        # --------------------------------------------------------------
        # Add additional constraints before the grouping, e.g. only MS2.
        # --------------------------------------------------------------
        _cnstr = []

        if only_with_rt:
            _cnstr.append("retention_time IS NOT NULL")
        if only_ms2:
            _cnstr.append("ms_level IS 'MS2'")
        if only_lc:
            _cnstr.append("instrument_type LIKE 'LC-%'")

        if len(_cnstr) > 0:
            _stmt.append("WHERE " + " AND ".join(_cnstr))

        # --------------------
        # Group the accessions
        # --------------------
        _stmt.append("GROUP BY ion_mode, instrument, instrument_type, column_name, column_temperature, flow_gradient, "
                     "  flow_rate, solvent_A, solvent_B, fragmentation_mode")

        cur = filter_db_conn.execute(" ".join(_stmt))

        # ============================================
        # Include the accessions as separate datasets:
        # ============================================
        for acc_pref_idx, row in enumerate(cur):
            dataset_identifier = "%s_%03d" % (accession_prefix, acc_pref_idx)

            for acc in sorted(row[1].split(",")):
                # -----------------------------------------------------------------------------
                # Update the molecule structure information in the spectrum file using PubChem:
                # -----------------------------------------------------------------------------
                if use_pubchem_structure_info and \
                        not specs[acc].update_molecule_structure_information_using_pubchem(pc_dbfn):
                    continue

                # -----------------------------------
                # Insert the spectrum to the database
                # -----------------------------------
                specs[acc].set("exact_mass_error_ppm",
                               get_mass_error_in_ppm(float(specs[acc].get("monoisotopic_mass")),
                                                     float(specs[acc].get("precursor_mz")),
                                                     specs[acc].get("precursor_type")))
                if specs[acc].get("exact_mass_error_ppm") is None:
                    LOGGER.info("{} Could not determine exact mass error. (precursor-type={})"
                                .format(acc, specs[acc].get("precursor_type")))
                    continue

                try:
                    with self._mb_conn:
                        self.insert_spectrum(dataset_identifier, contributor, specs[acc])
                except sqlite3.IntegrityError as err:
                    LOGGER.error("({}) SQLite integrity error: {}".format(acc, err))
                    continue

            # ----------------------------------------------------
            # Determine the number of spectra and unique compounds
            # ----------------------------------------------------
            with self._mb_conn:
                self._update_num_spectra_and_compounds_in_dataset_table(dataset_identifier)

        filter_db_conn.close()

    def insert_spectrum(self, dataset_identifier: str, contributor: str, spectrum):
        """

        :param spectrum:
        :return:
        """
        # ===========================
        # Insert Molecule Information
        # ===========================
        self._mb_conn.execute("INSERT OR IGNORE INTO molecules VALUES (%s)" % self._get_db_value_placeholders(11),
                              (
                                   spectrum.get("pubchem_id"),
                                   spectrum.get("inchi"),
                                   spectrum.get("inchikey"),
                                   spectrum.get("inchikey").split("-")[0],
                                   spectrum.get("inchikey").split("-")[1],
                                   spectrum.get("smiles_iso"),
                                   spectrum.get("smiles_can"),
                                   spectrum.get("exact_mass"),
                                   spectrum.get("monoisotopic_mass"),
                                   spectrum.get("molecular_formula"),
                                   spectrum.get("xlogp3")
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
        self._mb_conn.execute("INSERT OR IGNORE INTO datasets VALUES (%s)" % self._get_db_value_placeholders(19),
                              (
                                   dataset_identifier,
                                   contributor,
                                   ion_mode,
                                   -1,
                                   -1,
                                   spectrum.get("copyright"),
                                   spectrum.get("license"),
                                   spectrum.get("column_name"),
                                   None,  # TODO: column type, e.g. RP or HILIC, is not specified in the Massbank file
                                   spectrum.get("column_temperature"),
                                   spectrum.get("flow_gradient"),
                                   spectrum.get("flow_rate"),
                                   spectrum.get("solvent_A"),
                                   spectrum.get("solvent_B"),
                                   spectrum.get("solvent"),
                                   spectrum.get("instrument_type"),
                                   spectrum.get("instrument"),
                                   spectrum.get("fragmentation_mode"),
                                   None  # TODO: the column dead-time needs to be estimated from the data
                               ))

        # ===============================
        # Insert Spectra Meta Information
        # ===============================
        # Note: We can get a "FOREIGN KEY Constrained" SQLite Error here, if the dataset wasn't inserted due to missing
        #       information (e.g.).
        self._mb_conn.execute("INSERT INTO spectra_meta VALUES (%s)" % self._get_db_value_placeholders(11),
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
                                   spectrum.get("exact_mass_error_ppm")
                               ))

        # ====================
        # Insert Spectra Peaks
        # ====================
        _mz, _int = spectrum.get_mz(), spectrum.get_int()
        _n_peaks = len(_mz)
        _acc = [spectrum.get("accession")] * _n_peaks
        self._mb_conn.executemany("INSERT INTO spectra_peaks VALUES(%s)" % self._get_db_value_placeholders(3),
                                  zip(_acc, _mz, _int))

        # =====================
        # Insert Retention Time
        # =====================
        if spectrum.get("retention_time_unit") == '':
            if spectrum.get("retention_time") > 100:
                LOGGER.warning("(%s) Retention time unit default (=min) does not look reasonable: rt=%.3f"
                               % (spectrum.get("accession"), spectrum.get("retention_time")))

            self._mb_conn.execute("INSERT INTO spectra_raw_rts (spectrum, retention_time) VALUES(?, ?)",
                                  (
                                    spectrum.get("accession"),
                                    spectrum.get("retention_time"),
                                   ))
        else:
            self._mb_conn.execute("INSERT INTO spectra_raw_rts VALUES(?, ?, ?)",
                                  (
                                    spectrum.get("accession"),
                                    spectrum.get("retention_time"),
                                    spectrum.get("retention_time_unit")
                                   ))

    def iter_spectra(self, dataset, grouped=True, return_candidates=False, ppm=5, pc_dbfn=None):
        if grouped:
            rows = self._mb_conn.execute("SELECT GROUP_CONCAT(accession), dataset, molecule, "
                                         "       GROUP_CONCAT(precursor_mz), precursor_type, "
                                         "       GROUP_CONCAT(collision_energy), ms_type, "
                                         "       GROUP_CONCAT(resolution), fragmentation_mode "
                                         "    FROM spectra_meta "
                                         "    WHERE dataset IS ? "
                                         "    GROUP BY molecule, precursor_type, fragmentation_mode, ms_type",
                                         (dataset, ))
        else:
            rows = self._mb_conn.execute("SELECT accession, dataset, molecule, precursor_mz, precursor_type, "
                                         "       collision_energy, ms_type, resolution, fragmentation_mode "
                                         "   FROM spectra_meta "
                                         "   WHERE dataset IS ?", (dataset, ))

        # Open a connection to a local PubChemDB if needed and provided
        if return_candidates in ["mf", "mz"]:
            if pc_dbfn is None:
                raise ValueError("No local PubChem DB. Candidates cannot be queried.")

            pc_conn = sqlite3.connect("file:" + pc_dbfn + "?mode=ro", uri=True)  # open read-only
        else:
            pc_conn = None

        for row in rows:
            specs = []

            # -------------------------
            # Load compound information
            # -------------------------
            # Note: This is equal for all spectra, if grouped, as the compound is a grouping criteria.
            mol = self._mb_conn.execute("SELECT * FROM molecules WHERE cid = ?", (row[2],)).fetchall()[0]

            # --------------------------
            # Load candidate information
            # --------------------------
            if return_candidates == "mf":
                cands = pd.read_sql("SELECT * FROM compounds WHERE molecular_formula IS '%s'" % mol[8], con=pc_conn)
            elif return_candidates == "mz":
                min_exact_mass, max_exact_mass = self._get_ppm_window(mol[7], ppm)
                cands = pd.read_sql("SELECT * FROM compounds WHERE exact_mass BETWEEN %f AND %f" %
                                    (min_exact_mass, max_exact_mass), con=pc_conn)
            elif not return_candidates:
                cands = None
            else:
                raise ValueError("Invalid")

            # -----------------------------------------------------
            # Create a Spectrum object for all spectra in the group
            # -----------------------------------------------------
            for acc, ce in zip(row[0].split(","), row[5].split(",")):
                specs.append(massbank2db.spectrum.MBSpectrum())

                # --------------------------
                # Data about the acquisition
                # --------------------------
                specs[-1].set("accession", acc)
                specs[-1].set("dataset", row[1])
                specs[-1].set("precursor_mz", row[3])
                specs[-1].set("precursor_type", row[4])
                specs[-1].set("collision_energy", ce)
                specs[-1].set("ms_type", row[6])
                specs[-1].set("resolution", row[7])
                specs[-1].set("fragmentation_mode", row[8])

                # -------------------
                # Retention time data
                # -------------------
                rt, rt_unit = self._mb_conn.execute("SELECT retention_time, retention_time_unit FROM spectra_raw_rts"
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
                peaks = self._mb_conn.execute("SELECT mz, intensity FROM spectra_peaks "
                                              "    WHERE spectrum IS ? ORDER BY mz", (acc, ))
                specs[-1]._mz, specs[-1]._int = [], []
                for peak in peaks:
                    specs[-1]._mz.append(peak[0])
                    specs[-1]._int.append(peak[1])

            if not grouped:
                specs = specs[0]

            yield mol, specs, cands

        if pc_conn:
            pc_conn.close()

    def filter_datasets_by_number_of_unique_compounds(self, min_number_of_unique_compounds_per_dataset: int = 50):
        """

        :param min_number_of_unique_compounds_per_dataset:
        :return:
        """
        for ds in self._get_list_of_datasets():
            _num_cmp = self._mb_conn.execute(
                "SELECT COUNT(distinct molecule) FROM spectra_meta WHERE dataset IS ?", (ds, )).fetchone()[0]

            if _num_cmp < min_number_of_unique_compounds_per_dataset:
                LOGGER.info("Remove {} due to low number of unique compounds: {} < {}."
                            .format(ds, _num_cmp, min_number_of_unique_compounds_per_dataset))

                with self._mb_conn:
                    self._mb_conn.execute("DELETE FROM datasets WHERE name IS ?", (ds, ))  # triggers an SQLite event

        return self

    def filter_spectra_by_mass_error(self, max_exact_mass_error_ppm: float = 20):
        """

        :param max_exact_mass_error_ppm:
        :return:
        """
        with self._mb_conn:
            self._mb_conn.execute("DELETE FROM spectra_meta WHERE exact_mass_error_ppm > ?", (max_exact_mass_error_ppm, ))

        self._update_num_spectra_and_compounds_in_dataset_table()

        return self

    def filter_spectra_by_retention_time(self):
        pass

    def _get_list_of_datasets(self) -> List[str]:
        return [row[0] for row in self._mb_conn.execute("SELECT name FROM datasets")]

    def _update_num_spectra_and_compounds_in_dataset_table(self, dataset_identifier: Optional[str] = None) -> None:
        """
        Function to update the number of spectra and unique compounds for the given datasets.

        :param dataset_identifier: string or None, name of the dataset for which the statistics should be updated. If
            None, all datasets are updated.
        """
        if dataset_identifier is None:
            for ds in self._get_list_of_datasets():
                self._update_num_spectra_and_compounds_in_dataset_table(ds)
        else:
            _num_spec, _num_cmp = self._mb_conn.execute("SELECT COUNT(accession), COUNT(distinct molecule) "
                                                        "   FROM spectra_meta "
                                                        "   WHERE dataset IS ?", (dataset_identifier,)).fetchone()
            self._mb_conn.execute("UPDATE datasets "
                                  "   SET num_spectra = ?, num_compounds = ?"
                                  "   WHERE name IS ?", (_num_spec, _num_cmp, dataset_identifier))

    @staticmethod
    def _cands_to_metfrag_format(cands):
        cands_out = cands[["exact_mass", "InChI", "cid", "InChIKey", "molecular_formula", "SMILES_ISO"]]
        cands_out.assign(InChIKey1=cands_out["InChIKey"].apply(lambda _r: _r.split("-")[0]))
        cands_out.assign(InChIKey2=cands_out["InChIKey"].apply(lambda _r: _r.split("-")[1]))
        cands_out = cands_out.rename({"cid": "Identifier",
                                      "molecular_formula": "MolecularFormula",
                                      "exact_mass": "MonoisotopicMass",
                                      "SMILES_ISO": "SMILES"},
                                     axis=1)
        return cands_out.to_csv(sep="|", index=False)

    @staticmethod
    def _get_temporal_database(file_pth=":memory:"):
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
                     "  retention_time      FLOAT,"
                     "  fragmentation_mode  VARCHAR)")

        return conn

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

    # Load list of datasets provided by MassBank
    massbank_dir = "/home/bach/Documents/doctoral/data/MassBank-data_bachi55"
    mbds = pd.read_csv(os.path.join(massbank_dir, "List_of_Contributors_Prefixes_and_Projects.md"),
                       sep="|", skiprows=2, header=None) \
        .iloc[:, [1, 4]] \
        .applymap(str.strip) \
        .rename({1: "Contributor", 4: "AccPref"}, axis=1)  # type: pd.DataFrame

    # Filename of the MassBank (output) DB
    mb_dbfn = "tests/test_DB.sqlite"
    with MassbankDB(mb_dbfn) as mbdb:
        mbdb.initialize_tables(reset=True)

    # Filename of the local PubChem DB
    pc_dbfn = "/home/bach/Documents/doctoral/projects/local_pubchem_db/db_files/pubchem_01-02-2021.sqlite"

    # Insert spectra into the MassBank DB
    for idx, row in mbds.iterrows():
        LOGGER.info("== Process Dataset: {} ({} / {}) ==".format(row["Contributor"], idx + 1, len(mbds)))

        for pref in map(str.strip, row["AccPref"].split(",")):
            LOGGER.info("-- Process prefix: {}".format(pref))
            with MassbankDB(mb_dbfn) as mbdb:
                mbdb.insert_dataset(pref, row["Contributor"], massbank_dir, pc_dbfn=pc_dbfn,
                                    use_pubchem_structure_info=True)

    # Filter the database
    with MassbankDB(mb_dbfn) as mbdb:
        mbdb.filter_spectra_by_mass_error(max_exact_mass_error_ppm=20) \
            .filter_datasets_by_number_of_unique_compounds(min_number_of_unique_compounds_per_dataset=50)
