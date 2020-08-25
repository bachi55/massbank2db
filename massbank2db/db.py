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
        conn.execute(
            "CREATE TABLE datasets( \
                name                VARCHAR PRIMARY KEY  NOT NULL, \
                contributor         VARCHAR NOT NULL, \
                retention_time_unit VARCHAR NOT NULL, \
                copyright           VARCHAR, \
                license             VARCHAR, \
                column_name         VARCHAR, \
                column_temperature  VARCHAR, \
                flow_gradient       VARCHAR, \
                flow_rate           VARCHAR, \
                solvent_A           VARCHAR, \
                solvent_B           VARCHAR, \
                solvent             VARCHAR, \
                instrument_type     VARCHAR, \
                instrument          VARCHAR, \
             PRIMARY KEY())"
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
                retention_time      FLOAT NOT NULL, \
                collision_energy    FLOAT, \
                ms_level            INTEGER, \
                resolution          FLOAT, \
                fragmentation_type  VARCHAR, \
                origin              VARCHAR, \
             FOREIGN KEY(molecule)  REFERENCES molecules(cid),\
             FOREIGN KEY(dataset)   REFERENCES datasets(name) ON DELETE CASCADE)"
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


class MassbankDB(object):
    """

    """
    def __init__(self, file_pth):
        self.__conn = sqlite3.connect(file_pth)

    def __enter__(self):
        return self

    def __exit__(self, ext_type, exc_value, traceback):
        if isinstance(exc_value, Exception):
            self.__conn.rollback()

        self.__conn.close()

    def insert_dataset(self):
        pass

if __name__ == "__main__":
    create_db("tests/test_DB.sqlite")
