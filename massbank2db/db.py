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

from massbank2db.spectrum import MBSpectrum


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
                 inchikey3         VARCHAR NOT NULL, \
                 smiles_iso        VARCHAR NOT NULL, \
                 smiles_can        VARCHAR NOT NULL, \
                 molecular_weight  FLOAT NOT NULL, \
                 exact_mass        FLOAT NOT NULL, \
                 molecular_formula VARCHAR NOT NULL)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey_index ON molecules(inchikey)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey1_index ON molecules(inchikey1)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_exact_mass_index ON molecules(exact_mass)")
        conn.execute("CREATE INDEX IF NOT EXISTS molecules_mf_index ON molecules(molecular_formula)")

        # Datasets Meta-information table
        #   primary key: Accession prefix + some running id
        #                E.g. AU_001, AU_002
        conn.execute(
            "CREATE TABLE datasets( \
                name                VARCHAR PRIMARY KEY NOT NULL, \
                contributor         VARCHAR NOT NULL, \
                copyright           VARCHAR, \
                license             VARCHAR, \
                column_name         VARCHAR, \
                column_type         VARCHAR NOT NULL, \
                column_temperature  FLOAT, \
                flow_gradient       VARCHAR, \
                flow_rate           FLOAT, \
                solvent_A           VARCHAR, \
                solvent_B           VARCHAR, \
                solvent             VARCHAR, \
                instrument_type     VARCHAR, \
                instrument          VARCHAR, \
                column_dead_time    FLOAT)"
        )

        # Spectra Meta-information table
        conn.execute(
            "CREATE TABLE spectra_meta( \
                accession           VARCHAR PRIMARY KEY NOT NULL, \
                dataset             VARCHAR NOT NULL, \
                record_title        VARCHAR NOT NULL, \
                molecule            INTEGER NOT NULL, \
                ion_mode            VARCHAR NOT NULL, \
                precursor_mz        FLOAT NOT NULL, \
                precursor_type      FLOAT NOT NULL, \
                collision_energy    FLOAT, \
                ms_type             VARCHAR, \
                resolution          FLOAT, \
                fragmentation_type  VARCHAR, \
                origin              VARCHAR, \
             FOREIGN KEY(molecule)  REFERENCES molecules(cid),\
             FOREIGN KEY(dataset)   REFERENCES datasets(name) ON DELETE CASCADE)"
        )

        # Spectra Retention Time information table
        conn.execute(
            "CREATE TABLE spectra_raw_rts( \
                spectrum            VARCHAR NOT NULL, \
                retention_time      FLOAT NOT NULL, \
                retention_time_unit VARCHAR DEFAULT 'min', \
             FOREIGN KEY(spectrum) REFERENCES spectra_meta(accession) ON DELETE CASCADE)"
        )

        # Spectra Peak Information table
        conn.execute(
            "CREATE TABLE spectra_peaks( \
                mz        FLOAT NOT NULL, \
                intensity FLOAT NOT NULL, \
                spectrum  VARCHAR NOT NULL, \
             FOREIGN KEY(spectrum) REFERENCES spectra_meta(accession) ON DELETE CASCADE)")

        # Spectra candidate information
        conn.execute(
            "CREATE TABLE spectra_candidates( \
                spectrum  VARCHAR NOT NULL, \
                candidate INTEGER NOT NULL, \
                mf_equal  INTEGER NOT NULL, \
                mz_diff   FLOAT NOT NULL, \
             FOREIGN KEY(spectrum)  REFERENCES spectra_meta(accession)  ON DELETE CASCADE, \
             FOREIGN KEY(candidate) REFERENCES molecules(cid)           ON DELETE CASCADE, \
             PRIMARY KEY(spectrum, candidate))")


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
    """

    """
    def __init__(self, file_pth, only_with_rt=True, only_ms2=True, min_number_of_unique_compounds_per_dataset=50,
                 pubchem_file_pth=None):
        self.__conn = sqlite3.connect(file_pth)

        if pubchem_file_pth:
            self.__pc_conn = sqlite3.connect(pubchem_file_pth)
        else:
            self.__pc_conn = None

        if not (only_with_rt and only_ms2):
            raise NotImplementedError("The current filtering only includes accessions with MS2 and RTs.")

        # Different filters to exclude entries from the database
        self.__only_with_rt = only_with_rt
        self.__only_ms2 = only_ms2
        self.__min_number_of_unique_compounds_per_dataset = min_number_of_unique_compounds_per_dataset

    def __enter__(self):
        return self

    def __exit__(self, ext_type, exc_value, traceback):
        # Rollback Massbank database in case of an error
        if isinstance(exc_value, Exception):
            self.__conn.rollback()

        # Close connection to the Massbank database
        self.__conn.close()

        # Close connection to the pubchem database
        if self.__pc_conn:
            self.__pc_conn.close()

    def insert_dataset(self, accession_prefix, contributor, base_path, pubchem_):
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

            filter_db_conn.execute("INSERT INTO information VALUES(%s)" % (",".join(["?"] * len(new_row)),),
                                   new_row)

        # ================================================================
        # Build up a temporary database to filter and separate accessions:
        # ================================================================
        # - group accessions by their MS and LC configuration --> each combination becomes a separate dataset
        # - remove groups of accessions, where we have less than a certain number
        #   (see 'min_number_of_unique_compounds_per_dataset')
        filter_db_conn = get_temporal_database(".tmp2.sqlite")  # in memory
        specs = {}
        with filter_db_conn:
            for msfn in glob.iglob(os.path.join(base_path, contributor, accession_prefix + "[0-9]*.txt")):
                acc = os.path.basename(msfn).split(".")[0]
                specs[acc] = MBSpectrum(msfn)
                assert specs[acc].get("accession") == acc

                if not specs[acc].get("inchikey"):
                    print("WARNING: No InChIKey for '%s'." % os.path.basename(msfn))
                    continue

                _add_spec_to_filter_db(specs[acc])

        # =====================
        # Group the accessions:
        # =====================
        cur = filter_db_conn.execute(
            "SELECT COUNT(DISTINCT inchikey) AS num_unique_molecules, GROUP_CONCAT(accession) AS accessions \
                FROM information \
                WHERE retention_time IS NOT NULL AND ms_level IS 'MS2' \
                GROUP BY ion_mode, instrument, instrument_type, column_name, column_temperature, flow_gradient, \
                         flow_rate, solvent_A, solvent_B")

        # ============================================
        # Include the accessions as separate datasets:
        # ============================================
        for row in cur.fetchall():
            if row[0] < self.__min_number_of_unique_compounds_per_dataset:
                continue

            for acc in row[1].split(","):
                # -----------------------------------------------------------------------------
                # Update the molecule structure information in the spectrum file using Pubchem:
                # -----------------------------------------------------------------------------
                specs[acc].update_molecule_structure_information_using_pubchem(self.__pc_conn)


                "spectra_candidates", "spectra_peaks", "spectra_raw_rts", "spectra_meta", "datasets",
                "molecules"


        print("bla")


if __name__ == "__main__":
    dbfn = "tests/test_DB.sqlite"

    create_db(dbfn)

    with MassbankDB(dbfn) as mbdb:
        mbdb.insert_dataset("AU", "Athens_Univ", "/run/media/bach/EVO500GB/data/MassBank")