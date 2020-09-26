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
import unittest

from massbank2db.parser import get_AC_regex, get_CH_regex, get_meta_regex


class TestRegularExpressions(unittest.TestCase):
    def test_AC_retention_time(self):
        rex = get_AC_regex()["retention_time"]

        rt_str_val = [
            ("AC$CHROMATOGRAPHY: RETENTION_TIME 23.46 min (in paper: 23.3 min)", ("23.46", "min")),  # Chubu_Univ/UT001973.txt
            ("AC$CHROMATOGRAPHY: RETENTION_TIME 2.75", ("2.75", "")),  # AAFC/AC000338.txt
            ("AC$CHROMATOGRAPHY: RETENTION_TIME 9.367 min  ", ("9.367", "min")),  # Athens_Univ/AU592019.txt
            ("AC$CHROMATOGRAPHY: RETENTION_TIME  478.2 sec", ("478.2", "sec")),  # BS/BS001044.txt:
            ("AC$CHROMATOGRAPHY: RETENTION_TIME 618.594  sec", ("618.594", "sec")),  # OUF00490.txt
            ("AC$CHROMATOGRAPHY: RETENTION_TIME 618.594  s", ("618.594", "s"))
        ]

        for line, ref in rt_str_val:
            matches = rex[0].match(line)

            if ref:
                rt, rt_unit = ref
                self.assertEqual(2, len(matches.groups()))
                self.assertEqual(rt, matches.groups()[0])
                self.assertEqual(rt_unit, matches.groups()[1])
            else:
                self.assertEqual(None, matches)

    def test_CH_pubchem_id(self):
        rex = get_CH_regex()["pubchem_id"]

        cid_str_val = [("CH$LINK: PUBCHEM CID:8655", "8655"),
                       ("CH$LINK: PUBCHEM 3437", None),
                       ("CH$LINK: PUBCHEM SID:12293 CID:5281672", "5281672"),
                       ("CH$LINK: PUBCHEM CID:5280459 SID:4883", "5280459"),
                       ("CH$LINK: PUBCHEM CID 53297381", None),
                       ("CH$LINK: PUBCHEM CID 53297 3821", None),
                       ("CH$LINK: PUBCHEM SID 53297 3821", None),
                       ("CH$LINK: PUBCHEM cid:5280459", None),
                       ("CH$LINK: PUBCHEM SID:5", None),
                       ("CH$LINK: PUBCHEM CID:5", "5")]

        for line, cid in cid_str_val:
            for _rex in rex:
                matches = _rex.match(line)

                # We stop searching, ones a pattern has matched
                if matches:
                    break

            if cid:
                self.assertEqual(1, len(matches.groups()))
                self.assertEqual(cid, matches.groups()[0])
            else:
                self.assertIs(None, matches)

    def test_meta_accession(self):
        rex = get_meta_regex()["accession"]

        acc_str_val = [("ACCESSION: XY000010", "XY000010"),
                       ("ACCESSION: AUD43241", "AUD43241")]

        for line, acc in acc_str_val:
            matches = rex[0].match(line)

            self.assertEqual(acc, matches.groups()[0])


if __name__ == '__main__':
    unittest.main()
