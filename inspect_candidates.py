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

from massbank2db.db import MassbankDB
from massbank2db.spectrum import MBSpectrum

if __name__ == "__main__":
    mb_dbfn = "./massbank2db/tests/test_DB.sqlite"
    pc_dbfn = "/run/media/bach/EVO500GB/data/pubchem_24-06-2019/db/pubchem.sqlite"

    with MassbankDB(mb_dbfn) as mb_db:
        for mol, specs, cands in mb_db.iter_spectra(dataset="ET_000", return_candidates=False):
            print([spec.get("accession") for spec in specs], mol[3])

            MBSpectrum.merge_spectra(specs)

            print("bla")