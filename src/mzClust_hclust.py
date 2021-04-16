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
import numpy as np

from ctypes import *
from scipy.spatial.distance import pdist

# Compile the c-file
# gcc -c -fPIC mzClust_hclust.c -o mzClust_hclust.o
# gcc -shared -Wl,-soname,libclust.so -o libclust.so mzClust_hclust.o

if __name__ == "__main__":
    lib = cdll.LoadLibrary("libclust.so")

    x = [32, 4, 4, 6, 9]
    # TODO: Check pdist outputs the right dimension
    # TODO: Why do we need to use euclidean distance --> check in the original publication.
    d = pdist(np.array(x)[:, np.newaxis])
    g = [-1] * len(x)

    C_x = (c_double * len(x))(*x)
    C_d = (c_double * len(d))(*d)
    C_n = c_int(len(x))
    C_g = (c_int * len(x))(*g)
    C_eppm = c_double(0.01)
    C_eabs = c_double(0.01)

    print(d)
    print(C_x)

    lib.R_mzClust_hclust(byref(C_x), byref(C_n), byref(C_d), byref(C_g), byref(C_eppm), byref(C_eabs))

    for i in range(len(x)):
        print(C_g[i])