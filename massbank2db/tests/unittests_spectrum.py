####
#
# The 'massbank2db' package can be used to an build SQLite DB from the Massbank MS/MS repository.
#
#     Copyright (C) 2020 - 2022  Eric Bach <eric.bach@aalto.fi>
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
import unittest
import glob
import os
import numpy as np
import pandas as pd

from massbank2db.spectrum import MBSpectrum


class TestMBSpectrumInfoSanitizer(unittest.TestCase):
    def test_rt_sanitizer(self):
        out = MBSpectrum._sanitize_meta_information({"retention_time": ("430", "")})
        self.assertEqual(430, out["retention_time"])
        self.assertEqual(None, out["retention_time_unit"])

        out = MBSpectrum._sanitize_meta_information({"retention_time": ("430.3", "min")})
        self.assertEqual(430.3, out["retention_time"])
        self.assertEqual("min", out["retention_time_unit"])

        out = MBSpectrum._sanitize_meta_information({"retention_time": ("32", "s")})
        self.assertEqual(32, out["retention_time"])
        self.assertEqual("sec", out["retention_time_unit"])

        out = MBSpectrum._sanitize_meta_information({"retention_time": ("32", "m")})
        self.assertEqual(32, out["retention_time"])
        self.assertEqual("min", out["retention_time_unit"])


class TestMBSpectrumParsing(unittest.TestCase):
    def test_peak_parsing(self):
        # Spectrum 1
        peaks = [(65.0388, 454422.7),
                 (105.0703, 3594629.8),
                 (115.0544, 916283.5),
                 (121.0287, 482029.8),
                 (139.0544, 904286),
                 (163.0549, 764499.2),
                 (164.0627, 2677649.4),
                 (165.0706, 293933491.2),
                 (166.0782, 235412982.7),
                 (183.0812, 3256920.2),
                 (193.0767, 1661789.6),
                 (199.0315, 5538869.4),
                 (201.0474, 7384944.1)]
        spec = MBSpectrum("./example_massbank_records/EQ308406.txt")
        self.assertEqual(peaks, spec.get_peaks())

        # Spectrum 2
        peaks = [
            (116.050500, 2998.000000),
            (117.054700, 236.000000),
            (118.029100, 170.000000),
            (131.037600, 1.23552e13),
            (132.045500, 1241.000000),
            (133.048300, 116.000000),
            (159.032400, 1732.000000),
            (160.040200, 10392.000000),
            (161.043200, 1399.000000)]
        spec = MBSpectrum("./example_massbank_records/FIO00665.txt")
        self.assertEqual(peaks, spec.get_peaks())


class TestMBSpectrumToToolFormat(unittest.TestCase):
    def test_to_sirius(self):
        # Spectrum 1 -------------------
        spec = MBSpectrum("./example_massbank_records/FIO00665.txt")
        out = spec.to_sirius_format()
        self.assertIn("FIO00665.ms", out)
        self.assertIn("FIO00665.tsv", out)
        self.assertIsNone(out["FIO00665.tsv"])
        self.assertIn(">profile qtof", out["FIO00665.ms"])

        # check fragmentation peaks
        tmp = out["FIO00665.ms"].split("\n")
        for idx, _peak in enumerate(spec.get_peaks(), start=tmp.index(">ms2merged") + 1):
            _mz, _int = tmp[idx].split(" ")
            self.assertEqual(_peak, (float(_mz), float(_int)))

        # Spectrum 2 -------------------
        spec = MBSpectrum("./example_massbank_records/EQ308406.txt")
        out = spec.to_sirius_format()
        self.assertIn("EQ308406.ms", out)
        self.assertIn("EQ308406.tsv", out)
        self.assertIsNone(out["EQ308406.tsv"])
        self.assertIn(">profile orbitrap", out["EQ308406.ms"])

        # check fragmentation peaks
        tmp = out["EQ308406.ms"].split("\n")
        for idx, _peak in enumerate(spec.get_peaks(), start=tmp.index(">ms2merged") + 1):
            _mz, _int = tmp[idx].split(" ")
            self.assertEqual(_peak, (float(_mz), float(_int)))

        # Spectrum 3 --------------------
        spectra = []
        acc = []
        spec_cnt = 0
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA0004[01][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))
            acc.append(spectra[-1].get("accession"))
            spec_cnt += 1

        self.assertIn(
            ">ms2merged",
            MBSpectrum.merge_spectra(spectra, merge_peak_lists=True).to_sirius_format()["EA33002987.ms"]
        )

        self.assertEqual(
            spec_cnt,
            MBSpectrum.merge_spectra(spectra, merge_peak_lists=False)
                .to_sirius_format()["EA33002987.ms"]
                .count(">ms2peaks")
        )

    def test_to_sirius__gt_molecular_formula(self):
        acc = "EQ308406"
        self.assertIn("#formula",
                      MBSpectrum("./example_massbank_records/%s.txt" % acc).to_sirius_format()[acc + ".ms"])
        self.assertNotIn(">formula",
                         MBSpectrum("./example_massbank_records/%s.txt" % acc).to_sirius_format()[acc + ".ms"])
        self.assertNotIn("#formula",
                         MBSpectrum("./example_massbank_records/%s.txt" % acc)
                            .to_sirius_format(add_gt_molecular_formula=True)[acc + ".ms"])
        self.assertIn(">formula",
                      MBSpectrum("./example_massbank_records/%s.txt" % acc)
                         .to_sirius_format(add_gt_molecular_formula=True)[acc + ".ms"])

    def test_to_sirius__custom_candidate_db(self):
        # Spectrum 1 -------------------
        spec = MBSpectrum("./example_massbank_records/FIO00665.txt")
        out = spec.to_sirius_format(
            molecular_candidates=pd.read_csv("./example_massbank_records/FIO00665.tsv", sep="\t"))
        self.assertIn("FIO00665.ms", out)
        self.assertIn("FIO00665.tsv", out)
        self.assertIsNotNone(out["FIO00665.tsv"])

    def test_to_metfrag(self):
        # Spectrum 1 --------------------
        out = MBSpectrum("./example_massbank_records/FIO00665.txt").to_metfrag_format(
            **{"MetFragScoreWeights": [0.8, 0.2],
               "MetFragScoreTypes": ["FragmenterScore", "PubChemNumberPatents"],
               "LocalDatabasePath": "/path/to/db",
               "ResultsPath": "/path/to/results",
               "NumberThreads": 4,
               "PeakListPath": "/path/to/peaks"}
        )

        self.assertIn("FIO00665.peaks", out)
        self.assertIn("FIO00665.conf", out)
        self.assertIn("PeakListPath=/path/to/peaks/FIO00665.peaks", out["FIO00665.conf"])
        self.assertIn("MetFragScoreWeights=0.8,0.2\n", out["FIO00665.conf"])
        self.assertIn("MetFragScoreTypes=FragmenterScore,PubChemNumberPatents\n", out["FIO00665.conf"])
        self.assertIn("PrecursorIonMode=-1\n", out["FIO00665.conf"])
        self.assertIn("IsPositiveIonMode=False\n", out["FIO00665.conf"])

        # Spectrum 2 --------------------
        out = MBSpectrum("./example_massbank_records/EQ308406.txt").to_metfrag_format(
            **{"MetFragScoreWeights": [1.0],
               "MetFragScoreTypes": ["FragmenterScore"],
               "LocalDatabasePath": "/path/to/db",
               "ResultsPath": "/path/to/results",
               "NumberThreads": 4,
               "PeakListPath": "/path/to/peaks"}
        )

        self.assertIn("EQ308406.peaks", out)
        self.assertIn("EQ308406.conf", out)
        self.assertIn("PeakListPath=/path/to/peaks/EQ308406.peaks", out["EQ308406.conf"])
        self.assertIn("MetFragScoreWeights=1.0\n", out["EQ308406.conf"])
        self.assertIn("MetFragScoreTypes=FragmenterScore\n", out["EQ308406.conf"])
        self.assertIn("PrecursorIonMode=1\n", out["EQ308406.conf"])
        self.assertIn("IsPositiveIonMode=True\n", out["EQ308406.conf"])

        # Spectrum 3 --------------------
        spectra = []
        acc = []
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA0004[01][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))
            acc.append(spectra[-1].get("accession"))

        merged_spectrum = MBSpectrum.merge_spectra(spectra)
        out = merged_spectrum.to_metfrag_format(
            **{"MetFragScoreWeights": [1.0],
               "MetFragScoreTypes": ["FragmenterScore"],
               "LocalDatabasePath": "/path/to/db",
               "ResultsPath": "/path/to/results",
               "NumberThreads": 4,
               "PeakListPath": "/path/to/peaks"}
        )

        peaks_fn = merged_spectrum.get("accession") + ".peaks"
        config_fn = merged_spectrum.get("accession") + ".conf"

        self.assertIn("PeakListPath=/path/to/peaks/" + peaks_fn, out[config_fn])
        self.assertIn("MetFragScoreWeights=1.0\n", out[config_fn])
        self.assertIn("MetFragScoreTypes=FragmenterScore\n", out[config_fn])
        self.assertIn("PrecursorIonMode=1\n", out[config_fn])
        self.assertIn("IsPositiveIonMode=True\n", out[config_fn])


class TestMBSpectrumMerging(unittest.TestCase):
    def test_metainformation_merging__FIO00665(self):
        # Apply merge function to a single spectrum
        mb_fn = os.path.join("example_massbank_records", "FIO00665.txt")
        spec = MBSpectrum(mb_fn)
        spectra = [spec]
        acc_ref = [os.path.basename(mb_fn).split(".")[0]]
        rt_ref = spec.get("retention_time")
        precmz_ref = spec.get("precursor_mz")
        recordtitle_ref = spec.get("record_title")
        ce_ref = spec.get("collision_energy")

        # -------------------
        # WITH RT AGGREGATION
        # -------------------
        merged_spectrum = MBSpectrum.merge_spectra(spectra, rt_agg_fun=np.min)  # type: MBSpectrum

        self.assertEqual("FBZONXHGGPHHIY-UHFFFAOYSA-N", merged_spectrum.get("inchikey"))
        self.assertEqual(acc_ref, merged_spectrum.get("original_accessions"))
        self.assertEqual("FIO", merged_spectrum.get("accession")[:3])
        self.assertEqual(MBSpectrum._get_new_accession_id(merged_spectrum.get("original_accessions")),
                         merged_spectrum.get("accession"))
        self.assertEqual(rt_ref, merged_spectrum.get("retention_time"))
        self.assertEqual(precmz_ref, merged_spectrum.get("precursor_mz"))
        self.assertEqual(recordtitle_ref, merged_spectrum.get("record_title"))
        self.assertEqual([ce_ref], merged_spectrum.get("collision_energy"))

    def test_metainformation_merging__EAX000401(self):
        # Load the list of spectra to merge: EA0004[01][0-9].txt --> EAX000401.txt
        spectra = []
        acc_ref = []
        rt_ref = []
        precmz_ref = []
        recordtitle_ref = []
        ce_ref = []
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA0004[01][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))

            # collect some reference meta-information
            acc_ref.append(os.path.basename(mb_fn).split(".")[0])
            rt_ref.append(spectra[-1].get("retention_time"))
            precmz_ref.append(spectra[-1].get("precursor_mz"))
            recordtitle_ref.append(spectra[-1].get("record_title"))
            ce_ref.append(spectra[-1].get("collision_energy"))

        # -------------------
        # WITH RT AGGREGATION
        # -------------------
        merged_spectrum = MBSpectrum.merge_spectra(spectra, rt_agg_fun=np.mean)  # type: MBSpectrum

        self.assertEqual("OUSYWCQYMPDAEO-UHFFFAOYSA-N", merged_spectrum.get("inchikey"))
        self.assertEqual(acc_ref, merged_spectrum.get("original_accessions"))
        self.assertEqual("EA", merged_spectrum.get("accession")[:2])
        self.assertEqual(MBSpectrum._get_new_accession_id(merged_spectrum.get("original_accessions")),
                         merged_spectrum.get("accession"))
        self.assertEqual("min", merged_spectrum.get("retention_time_unit"))
        self.assertEqual(np.mean(rt_ref), merged_spectrum.get("retention_time"))
        self.assertEqual(precmz_ref[0], merged_spectrum.get("precursor_mz"))
        self.assertEqual(recordtitle_ref, merged_spectrum.get("record_title"))
        self.assertEqual(ce_ref, merged_spectrum.get("collision_energy"))

        # ----------------------
        # WITHOUT RT AGGREGATION
        # ----------------------
        merged_spectrum = MBSpectrum.merge_spectra(spectra, rt_agg_fun=None)  # type: MBSpectrum

        self.assertEqual("min", merged_spectrum.get("retention_time_unit"))
        self.assertEqual(rt_ref, merged_spectrum.get("retention_time"))

    def test_spectra_merging__EAX000401(self):
        """
        We compare our spectra merging with the strategy applied in [1] and originally proposed in [2].

        References:
            [1] "MetFrag relaunched: incorporating strategies beyond in silico fragmentation" by Ruttkies et al. (2016)
            [2] "Alignment of high resolution mass spectra: development of a heuristic approach for metabolomics" by
                Kazmi et al. (2006)
        """
        # Load the list of spectra to merge: EA0004[01][0-9].txt --> EAX000401.txt
        spectra = []
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA0004[01][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))

        # Run the spectra merging using hierarchical clustering
        merged_spectrum = MBSpectrum.merge_spectra(spectra)  # type: MBSpectrum

        # Merged spectrum as used by [1]
        peaks_ref = [
            (53.03852, 55),
            (57.0447285714286, 93),
            (65.0386, 116),
            (77.0386, 884),
            (81.03345, 12),
            (85.0396545454545, 17),
            (91.0542714285714, 123),
            (92.0494875, 377),
            (95.04925, 112),
            (103.041733333333, 15),
            (104.049541666667, 999),
            (105.044757142857, 271),
            (105.069975, 12),
            (110.060033333333, 7),
            (119.060475, 631),
            (130.04005, 11),
            (130.0652, 49),
            (131.07295, 24),
            (142.0652, 9),
            (147.0554, 3),
            (160.087069230769, 999),
            (188.082038461538, 999)
        ]

        mzs_ref = list(zip(*peaks_ref))[0]
        ints_ref = np.array(list(zip(*peaks_ref))[1])

        np.testing.assert_almost_equal(mzs_ref, merged_spectrum.get_mz())
        np.testing.assert_almost_equal(ints_ref / 999, merged_spectrum.get_int(), decimal=3)

    def test_spectra_merging__no_normalization__EAX000401(self):
        """
        We compare our spectra merging with the strategy applied in [1] and originally proposed in [2].

        References:
            [1] "MetFrag relaunched: incorporating strategies beyond in silico fragmentation" by Ruttkies et al. (2016)
            [2] "Alignment of high resolution mass spectra: development of a heuristic approach for metabolomics" by
                Kazmi et al. (2006)
        """
        # Load the list of spectra to merge: EA0004[01][0-9].txt --> EAX000401.txt
        spectra = []
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA0004[01][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))

        # Run the spectra merging using hierarchical clustering
        merged_spectrum = MBSpectrum.merge_spectra(spectra, normalize_peaks_before_merge=False)  # type: MBSpectrum

        # Merged spectrum as used by [1]
        peaks_ref = [
            (53.03852, 55),
            (57.0447285714286, 93),
            (65.0386, 116),
            (77.0386, 884),
            (81.03345, 12),
            (85.0396545454545, 17),
            (91.0542714285714, 123),
            (92.0494875, 377),
            (95.04925, 112),
            (103.041733333333, 15),
            (104.049541666667, 999),
            (105.044757142857, 271),
            (105.069975, 12),
            (110.060033333333, 7),
            (119.060475, 631),
            (130.04005, 11),
            (130.0652, 49),
            (131.07295, 24),
            (142.0652, 9),
            (147.0554, 3),
            (160.087069230769, 999),
            (188.082038461538, 999)
        ]

        mzs_ref = list(zip(*peaks_ref))[0]

        np.testing.assert_almost_equal(mzs_ref, merged_spectrum.get_mz())

    def test_spectra_merging__EAX281502(self):
        """
        We compare our spectra merging with the strategy applied in [1] and originally proposed in [2].

        References:
            [1] "MetFrag relaunched: incorporating strategies beyond in silico fragmentation" by Ruttkies et al. (2016)
            [2] "Alignment of high resolution mass spectra: development of a heuristic approach for metabolomics" by
                Kazmi et al. (2006)
        """
        # Load the list of spectra to merge: EA2815[56][0-9].txt --> EAX281502.txt
        spectra = []
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA2815[56][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))

        # Run the spectra merging using hierarchical clustering
        merged_spectrum = MBSpectrum.merge_spectra(spectra)  # type: MBSpectrum

        # Merged spectrum as used by [1]
        peaks_ref = [
            (81.0220667, 257.0000000),
            (96.0095000, 226.0000000),
            (109.0170500, 259.0000000),
            (111.0198667, 290.0000000),
            (116.0504625, 559.0000000),
            (118.0663111, 999.0000000),
            (126.1288286, 348.0000000),
            (156.0821000, 84.0000000),
            (170.1188000, 2.0000000),
            (172.0770500, 54.0000000),
            (174.0561667, 114.0000000),
            (182.1189250, 95.0000000),
            (197.1296571, 173.0000000),
            (200.0717900, 999.0000000),
            (227.1402000, 67.0000000),
            (230.1551000, 290.0000000),
            (244.0616500, 15.0000000),
            (257.2023400, 218.0000000),
            (273.1971500, 176.0000000),
            (276.0872500, 3.0000000),
            (283.1819500, 70.0000000),
            (285.1611800, 381.0000000),
            (301.1921875, 459.0000000),
            (317.1871200, 142.0000000),
            (327.1710000, 128.0000000),
            (331.2034500, 3.0000000),
            (345.1821167, 999.0000000),
            (359.1974333, 177.0000000),
            (377.2083750, 999.0000000)
        ]

        mzs_ref = list(zip(*peaks_ref))[0]
        ints_ref = np.array(list(zip(*peaks_ref))[1])

        np.testing.assert_almost_equal(mzs_ref, merged_spectrum.get_mz())
        np.testing.assert_almost_equal(ints_ref / 999, merged_spectrum.get_int(), decimal=3)

    def test_spectra_merging__EAX000402(self):
        """
        We compare our spectra merging with the strategy applied in [1] and originally proposed in [2].

        References:
            [1] "MetFrag relaunched: incorporating strategies beyond in silico fragmentation" by Ruttkies et al. (2016)
            [2] "Alignment of high resolution mass spectra: development of a heuristic approach for metabolomics" by
                Kazmi et al. (2006)
        """
        # self.skipTest("Intensities are wired in the original file: EAX000402")

        # Load the list of spectra to merge: EA0004[56][0-9].txt --> EAX000402.txt
        spectra = []
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA0004[56][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))

        # Run the spectra merging using hierarchical clustering
        merged_spectrum = MBSpectrum.merge_spectra(spectra)  # type: MBSpectrum

        # Merged spectrum as used by [1]
        peaks_ref = [
            (117.0347000, 999),
            (186.0677308, 999)
        ]

        mzs_ref = list(zip(*peaks_ref))[0]
        ints_ref = np.array(list(zip(*peaks_ref))[1])

        np.testing.assert_almost_equal(mzs_ref, merged_spectrum.get_mz())
        np.testing.assert_almost_equal(ints_ref / 999, merged_spectrum.get_int())

    def test_cfmid_output_parser(self):
        # Reference peaks for the example spectrum 2931.txt
        peaks_ref = [
            [
                (243.07, 6.2),
                (271.06010, 100.0)
            ],
            [
                (161.02332, 48.5),
                (163.03897, 24.9),
                (215.07027, 17.3),
                (243.06519, 25.8),
                (253.04954, 15.0),
                (271.06010, 100.0)
            ],
            [
                (65.03858, 39.7),
                (75.02293, 17.6),
                (77.03858, 34.4),
                (93.03349, 20.4),
                (109.02841, 25.7),
                (111.04406, 49.2),
                (113.05971, 19.1),
                (121.02841, 46.3),
                (123.04406, 23.2),
                (131.01276, 34.6),
                (133.02841, 33.4),
                (135.04406, 22.4),
                (137.02332, 100.0),
                (161.02332, 65.0),
                (163.03897, 26.2),
                (177.01824, 25.9),
                (187.03897, 21.4),
                (189.05462, 23.3),
                (197.02332, 18.0),
                (201.05462, 25.9),
                (211.03897, 38.6),
                (213.05462, 82.4),
                (215.07027, 30.6),
                (225.05462, 71.4),
                (241.04954, 15.4),
                (243.06519, 17.8),
                (253.04954, 17.9)
            ]
        ]

        # --- Load the insilico spectrum (each energy separately)
        spec = MBSpectrum.from_cfmid_output("example_cfmid_outputs/2931.txt", cfmid_4_format=True, merge_energies=False)

        for i in range(3):
            self.assertEqual("ID2931%d" % i, spec[i].get("accession"))
            self.assertEqual("energy%d" % i, spec[i].get("collision_energy"))

        for i in range(3):
            self.assertEqual(len(peaks_ref[i]), len(spec[i].get_peaks()))
            self.assertListEqual(peaks_ref[i], spec[i].get_peaks())

        # --- Load the insilico spectrum and merge the energies into a single spectrum
        spec = MBSpectrum.from_cfmid_output("example_cfmid_outputs/2931.txt", cfmid_4_format=True, merge_energies=True)

        self.assertIsInstance(spec, MBSpectrum)
        self.assertListEqual(["ID2931%d" % i for i in range(3)], spec.get("original_accessions"))

        # Peak only appears in one energy
        self.assertIn((65.03858, 39.7 / 100), spec.get_peaks())

        # Peak that appears in multiple energies
        self.assertIn((161.02332, 65.0 / 100), spec.get_peaks())


if __name__ == '__main__':
    unittest.main()
