####
#
# The 'massbank2db' package can be used to an build SQLite DB from the Massbank MS/MS repository.
#
#     Copyright (C) 2022  Eric Bach <eric.bach@aalto.fi>
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
import os
import logging
import pandas as pd
import argparse

from massbank2db.version import __version__ as massbank2db_version
from massbank2db.db import MassbankDB
from massbank2db.db import LOGGER as MB_LOGGER


if __name__ == "__main__":
    MB_LOGGER.setLevel(logging.INFO)

    # Set up and parse the input arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "massbank_repo_dir",
        type=str,
        help="Directory containing the local MassBank copy. It is assumed to be a clone of the MassBank GitHub "
             "repository."
    )
    arg_parser.add_argument("massbank_db_fn", type=str, help="Output filename of the MassBank SQLite database.")
    arg_parser.add_argument(
        "--use_pubchem_structure_info",
        action="store_true",
        help="Indicating whether a local PubChem copy should be used to update the molecular structure information."
    )
    arg_parser.add_argument(
        "--pubchem_db_fn",
        type=str,
        help="Filename of the local PubChem SQLite database file.",
        default=None
    )
    args = arg_parser.parse_args()

    # Load list of datasets provided by MassBank
    mbds = pd.read_csv(
        os.path.join(args.massbank_repo_dir, "List_of_Contributors_Prefixes_and_Projects.md"), sep="|", skiprows=2,
        header=None
    )
    mbds = mbds.iloc[:, [1, 4]].applymap(str.strip).rename({1: "Contributor", 4: "AccPref"}, axis=1)

    # Filename of the MassBank (output) DB
    _odir = os.path.dirname(args.massbank_db_fn)
    _ofn, _oext = os.path.basename(args.massbank_db_fn).split(os.extsep)
    _ofn = os.extsep.join([_ofn + "__%s" % massbank2db_version, _oext])
    massbank_db_fn = os.path.join(_odir, _ofn)
    MB_LOGGER.info("Output DB filename: %s" % massbank_db_fn)

    if os.path.exists(massbank_db_fn):
        raise RuntimeError("Output database file already exists: '%s'." % massbank_db_fn)

    with MassbankDB(massbank_db_fn) as mbdb:
        mbdb.initialize_tables(reset=True)

    # Insert spectra into the MassBank DB
    for idx, (_, row) in enumerate(mbds.iterrows()):
        MB_LOGGER.info("== Process Dataset: {} ({} / {}) ==".format(row["Contributor"], idx + 1, len(mbds)))

        for pref in map(str.strip, row["AccPref"].split(",")):
            MB_LOGGER.info("-- Process prefix: {}".format(pref))

            with MassbankDB(massbank_db_fn) as mbdb:
                mbdb.insert_dataset(
                    pref, row["Contributor"], args.massbank_repo_dir,
                    use_pubchem_structure_info=args.use_pubchem_structure_info, pc_dbfn=args.pubchem_db_fn,
                )

    # Filter the database
    with MassbankDB(massbank_db_fn) as mbdb:
        mbdb.filter_spectra_by_mass_error(max_exact_mass_error_ppm=20) \
            .filter_datasets_by_number_of_unique_compounds(min_number_of_unique_compounds_per_dataset=50)
