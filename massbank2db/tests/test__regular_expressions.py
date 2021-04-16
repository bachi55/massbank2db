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
import re


rt_regex = '^AC\$CHROMATOGRAPHY:\s+RETENTION.*TIME\s+(\d*[.,]?\d*)\s*($|min|sec)'

rt_strs = [
    "AC$CHROMATOGRAPHY: RETENTION_TIME 23.46 min (in paper: 23.3 min)",  # Chubu_Univ/UT001973.txt
    "AC$CHROMATOGRAPHY: RETENTION_TIME 2.75",                            # AAFC/AC000338.txt
    "AC$CHROMATOGRAPHY: RETENTION_TIME 9.367 min",                       # Athens_Univ/AU592019.txt
    "AC$CHROMATOGRAPHY: RETENTION_TIME 478.2 sec",                       # BS/BS001044.txt:
    "AC$CHROMATOGRAPHY: RETENTION_TIME 618.594  sec"                     # OUF00490.txt
]

for rt_str in rt_strs:
    print(rt_str)

    m = re.search(rt_regex, rt_str, re.IGNORECASE)
    print("\tMatch:", m)

    for i in [1, 2]:
        print("\t", m.group(i))

    print(m.groups())
