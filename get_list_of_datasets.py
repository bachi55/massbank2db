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

import pandas as pd
import argparse
import os
import glob
import sqlite3

from typing import Optional

from massbank2db.parser import parse_info, get_AC_regex


def create_temporal_database(file_pth=":memory:") -> sqlite3.Connection:
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
                 "  column_temperature  VARCHAR)")
    return conn


def _maintain_parser_output(d: dict, keys_to_keep: Optional[list] = None) -> tuple:
    out = tuple()

    if keys_to_keep:
        keys = keys_to_keep
    else:
        keys = d.keys()

    for k in keys:
        v = d[k]

        if not v:
            out += (None,)
        else:
            out += v

    return out


if __name__ == "__main__":
    # Handle input arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument("massbank_dir", help="Directory containing the MassBank data as provided in the Git "
                                                "repository: 'https://github.com/MassBank/MassBank-data'.")
    args = argparser.parse_args()

    # Load list of datasets provided by Massbank
    mbds = pd.read_csv(os.path.join(args.massbank_dir, "List_of_Contributors_Prefixes_and_Projects.md"),
                       sep="|", skiprows=2, header=None) \
        .iloc[:, [1, 4]] \
        .applymap(str.strip) \
        .rename({1: "Dataset", 4: "AccPref"}, axis=1)  # type: pd.DataFrame

    # Load all accession into a temporal database
    db_conn = create_temporal_database(".tmp.sqlite")

    for idx, row in mbds.iterrows():
        ds_dir = os.path.join(args.massbank_dir, row["Dataset"])

        for pref in row["AccPref"].split(","):
            pref = pref.strip()
            print(os.path.join(ds_dir, pref + "*.txt"))

            for fn in sorted(glob.glob(os.path.join(ds_dir, pref + "*.txt"))):
                info = parse_info(fn, regex=get_AC_regex())
                if info["solvent_A"] and info["solvent_B"]:
                    info["solvent"] = None

                info = _maintain_parser_output(info, keys_to_keep=["instrument_type", "instrument", "ion_mode",
                                                                   "column_name", "flow_gradient", "flow_rate",
                                                                   "solvent_A", "solvent_B", "solvent",
                                                                   "column_temperature"])

                with db_conn:
                    db_conn.execute("INSERT INTO information VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                    (fn.split("/")[-1].split(".")[0], row["Dataset"], pref) + info)

    db_conn.close()
