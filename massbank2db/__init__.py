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
from massbank2db.version import __version__

from os.path import dirname as __os_dirname__
from os.path import join as __os_join__
from ctypes import cdll
from platform import python_version_tuple

__python__version__ = python_version_tuple()
HCLUST_LIB = cdll.LoadLibrary(
    __os_join__(__os_dirname__(__file__),
                "mzClust_hclust.cpython-%s%s-x86_64-linux-gnu.so" % (__python__version__[0], __python__version__[1])))

