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

from massbank2db.db import MassbankDB
from massbank2db.spectrum import MBSpectrum

if __name__ == "__main__":
    mb_dbfn = "./massbank2db/tests/test_DB.sqlite"
    pc_dbfn = "/run/media/bach/EVO500GB/data/pubchem_24-06-2019/db/pubchem.sqlite"

    with MassbankDB(mb_dbfn) as mb_db:
        for mol, specs, cands in mb_db.iter_spectra(dataset="ET_000", return_candidates=False):
            print([spec.get("accession") for spec in specs], mol[3])

            # MBSpectrum.merge_spectra(specs)

            print(specs[0].get("precursor_mz"))

            print("bla")