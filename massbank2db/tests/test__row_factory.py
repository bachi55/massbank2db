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
import sqlite3


from massbank2db.utils import named_row_factory


if __name__ == "__main__":
    mb_dbfn = "test_DB.sqlite"

    con = sqlite3.connect(mb_dbfn)
    con.row_factory = named_row_factory

    res = con.execute("SELECT * FROM spectra_meta LIMIT 5")

    for row in res:
        print(row)

    con.row_factory = None
    res = con.execute("SELECT * FROM spectra_meta LIMIT 5")

    for row in res:
        print(row)
