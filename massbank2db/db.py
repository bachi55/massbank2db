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
                 "  retention_time      VARCHAR)")

    return conn


class MassbankDB(object):
    """

    """
    def __init__(self, file_pth, only_with_rt=True, only_ms2=False, min_number_of_unique_compounds_per_dataset=50):
        self.__conn = sqlite3.connect(file_pth)

        # Different filters to exclude entries from the database
        self.__only_with_rt = only_with_rt
        self.__only_ms2 = only_ms2
        self.__min_number_of_unique_compounds_per_dataset = min_number_of_unique_compounds_per_dataset

    def __enter__(self):
        return self

    def __exit__(self, ext_type, exc_value, traceback):
        if isinstance(exc_value, Exception):
            self.__conn.rollback()

        self.__conn.close()

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
        def _add_info_to_filter_db(filter_db_conn, info):
            acc = os.path.basename(msfn).split(".")[0]

            # Specify the solvent information
            solvent_A = info["solvent_A"]
            solvent_B = info["solvent_B"]
            if solvent_A and solvent_B:
                solvent = None
            else:
                solvent = info["solvent"]

            # Sanitize the RT information
            rt = info["retention_time"]
            if rt:
                rt = rt[0]

            new_row = (acc, contributor, accession_prefix) + info["instrument_type"] + info["instrument"] + \
                      info["ion_mode"] + info["column_name"] + info["flow_gradient"] + info["flow_rate"] + \
                      (solvent_A, solvent_B, solvent) + info["column_temperature"] + info["inchikey"] + \
                      info["ms_type"] + (rt, )

            filter_db_conn.execute("INSERT INTO information VALUES(%s)" % tuple(",".join(['?'] * len(new_row))),
                                   new_row)

        filter_db_conn = get_temporal_database()  # in memory

        for msfn in glob.iglob(os.path.join(base_path, contributor, accession_prefix + "[0-9]*.txt")):
            spec = MBSpectrum(msfn)

            if not spec._meta_information["inchikey"]:
                print("WARNING: No InChIKey for '%s'." % os.path.basename(msfn))
                continue

            _add_info_to_filter_db(filter_db_conn, spec)

            print("bla")


if __name__ == "__main__":
    dbfn = "tests/test_DB.sqlite"

    create_db(dbfn)

    with MassbankDB(dbfn) as mbdb:
        mbdb.insert_dataset("AU", "Athens_Univ", "/run/media/bach/EVO500GB/data/MassBank")