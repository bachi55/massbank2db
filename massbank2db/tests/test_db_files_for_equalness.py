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
