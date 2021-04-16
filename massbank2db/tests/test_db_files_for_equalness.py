####
#
# The 'massbank2db' package can be used to an build SQLite DB from the Massbank MS/MS repository.
#
#     Copyright (C) 2020 - 2021  Eric Bach <eric.bach@aalto.fi>
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
####
import pandas as pd
import sqlite3
import argparse
import sys


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("db1", type=str)
    arg_parser.add_argument("db2", type=str)
    arg_parser.add_argument("--tables_to_compare", nargs="+", default=["molecules", "datasets", "spectra_meta",
                                                                       "spectra_peaks", "spectra_raw_rts"])
    args = arg_parser.parse_args()

    # Open DB connections
    con1 = sqlite3.connect("file:" + args.db1 + "?mode=ro", uri=True)
    con2 = sqlite3.connect("file:" + args.db2 + "?mode=ro", uri=True)

    try:
        for tab in args.tables_to_compare:
            data1 = pd.read_sql_query("SELECT * FROM %s" % tab, con1)
            data2 = pd.read_sql_query("SELECT * FROM %s" % tab, con2)

            if not data1.equals(data2):
                raise ValueError("Not all data in table '%s' is equal." % tab)
    except ValueError as err:
        print(err)
        sys.exit(1)
    finally:
        # Close DB connections
        con1.close()
        con2.close()

    sys.exit(0)
