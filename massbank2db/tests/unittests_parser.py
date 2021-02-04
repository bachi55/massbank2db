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
import numpy as np

from massbank2db.parser import get_AC_regex, get_CH_regex, get_meta_regex
from massbank2db.parser import parse_column_name, parse_flow_rate_string


class TestParsingOfColumnInformation(unittest.TestCase):
    def test_column_dimensions_from_column_name(self):
        column_names = [
            "Agilent RRHD Eclipse 50 x 2 mm, 1.8 uM",
            "Acclaim RSLC C18 2.2um, 2.1x100mm, Thermo",
            "ACQUITY UPLC BEH Amide 1.7 um 2.1x100mm, Waters",
            "Acclaim RSLC C18 2.2um, 2.1x100mm, Thermo",
            "BEH C18 1.7um, 2.1x100mm, Waters",
            "Waters Acquity BEH C18 1.7um x 2.1 x 150 mm",
            "Kinetex C18 EVO 2.6 um, 2.1x50 mm, precolumn 2.1x5 mm, Phenomenex",
            "Develosil C30, Nomura Chemical",
            "XBridge C18 3.5um, 2.1x50mm, Waters",
            "Atlantis T3 3um, 3.0x150mm, Waters with guard column",
            "XBridge C18 3.5um, 2.1x150mm, Waters",
            "Acquity BEH C18 1.7um, 2.1x150mm (Waters)",
            "Symmetry C18 Column, Waters",
            "Kinetex Evo C18 2.6 um 50x2.1 mm, Phenomenex",
            "Acquity bridged ethyl hybrid C18 (1.7 um, 2.1 mm * 100 mm, Waters)",
            "Acquity UPLC Peptide BEH C18 column (50*2.1 mm; 1.7 um; 130A)(Waters Co.,Milford, MA, USA)",
            "Kinetex Core-Shell C18 2.6 um, 3.0 x 100 mm, Phenomenex",
            "Agilent C8 Cartridge Column 2.1X30mm 3.5 micron (guard); Agilent SB-Aq 2.1x50mm 1.8 micron (analytical)"
        ]
        diameters = [(2, "mm"), (2.1, ""), (2.1, ""), (2.1, ""), (2.1, ""), (2.1, ""), (2.1, ""), None, (2.1, ""),
                     (3, ""), (2.1, ""), (2.1, ""), None, (2.1, "mm"), (2.1, "mm"), (2.1, "mm"), (3.0, ""), (2.1, "")]
        lengths = [(50, ""), (100, "mm"), (100, "mm"), (100, "mm"), (100, "mm"), (150, "mm"), (50, "mm"), None,
                   (50, "mm"), (150, "mm"), (150, "mm"), (150, "mm"), None, (50, ""), (100, "mm"), (50, ""), (100, "mm"),
                   (50, "mm")]

        assert len(column_names) == len(diameters)
        assert len(column_names) == len(lengths)
        assert len(diameters) == len(lengths)

        for idx, cn in enumerate(column_names):
            dia, length = parse_column_name(cn)

            self.assertEqual(diameters[idx], dia)
            self.assertEqual(lengths[idx], length)

    def test_flow_rate_parsing(self):
        flow_rate_strs = [
            "0.3 mL min-1",
            "200 uL/min at 0-3 min, 400 uL/min at 14 min, 480 uL/min at 16-19 min, 200 uL/min at 19.1-20 min",
            "560 uL / min",
            "add 100 uL/min",
            "200 ul/min",
            "360 uL/min",
            "0.20 mL/min",
            "0.3 ml/min",
            "0.3 ml/min at 0 min, 0.3 ml/min at 10 min, 0.4 ml/min at 10.1 min, 0.4 ml/min at 14.4 min, 0.3 ml/min at 14.5 min",
            "0.20 mmL/min",
            "200 mL/min"]
        flow_rates = [(0.3, "mL/min"),
                      ([200, 400, 480, 200], "uL/min"),
                      (560, "uL/min"),
                      (100, "uL/min"),
                      (200, "ul/min"),
                      (360, "uL/min"),
                      (0.2, "mL/min"),
                      (0.3, "ml/min"),
                      ([0.3, 0.3, 0.4, 0.4, 0.3], "ml/min"),
                      None,
                      None]

        assert len(flow_rates) == len(flow_rate_strs)

        for idx, frstr in enumerate(flow_rate_strs):
            self.assertEqual(flow_rates[idx], parse_flow_rate_string(frstr))

    def test_flow_rate_parsing__aggregation(self):
        flow_rate_strs = [
            "200 uL/min at 0-3 min, 400 uL/min at 14 min, 480 uL/min at 16-19 min, 200 uL/min at 19.1-20 min",
            "0.3 ml/min",
            "0.3 ml/min at 0 min, 0.3 ml/min at 10 min, 0.4 ml/min at 10.1 min, 0.4 ml/min at 14.4 min, 0.3 ml/min at 14.5 min"]
        flow_rates = [([200, 400, 480, 200], "uL/min"),
                      (0.3, "ml/min"),
                      ([0.3, 0.3, 0.4, 0.4, 0.3], "ml/min")]

        assert len(flow_rates) == len(flow_rate_strs)

        for idx, frstr in enumerate(flow_rate_strs):
            # MEAN
            _val, _unit = parse_flow_rate_string(frstr, aggregate_flow_rates=np.mean)
            self.assertEqual(np.mean(flow_rates[idx][0]), _val)
            self.assertTrue(np.isscalar(_val))
            self.assertEqual(flow_rates[idx][1], _unit)

            # MIN
            _val, _unit = parse_flow_rate_string(frstr, aggregate_flow_rates=np.min)
            self.assertEqual(np.min(flow_rates[idx][0]), _val)
            self.assertTrue(np.isscalar(_val))
            self.assertEqual(flow_rates[idx][1], _unit)

    def test_flow_rate_parsing__TODO(self):
        self.skipTest("TODO: Cannot handle different units when multiple flow rates are reported.")

        flow_rate_strs = ["0.3 ml/min at 0 min, 200 ul/min"]
        flow_rates = [None]

        assert len(flow_rates) == len(flow_rate_strs)

        for idx, frstr in enumerate(flow_rate_strs):
            self.assertEqual(flow_rates[idx], parse_flow_rate_string(frstr))

    def test_flow_rate_parsing__currently_not_handled(self):
        self.skipTest("TODO")

        "200-320 (0-1 min); 200 (1-38 min) ul/min'"


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
