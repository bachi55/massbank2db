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
import glob
import os
import numpy as np

from massbank2db.spectrum import MBSpectrum


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
        self.assertEqual(peaks, spec.get_peak_list_as_tuples())

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
        self.assertEqual(peaks, spec.get_peak_list_as_tuples())


class TestMBSpectrumToToolFormat(unittest.TestCase):
    def test_to_metfrag(self):
        # Spectrum 1 --------------------
        out = MBSpectrum("./example_massbank_records/FIO00665.txt")._to_metfrag_format(
            **{"MetFragScoreWeights": [0.8, 0.2],
               "MetFragScoreTypes": ["FragmenterScore", "PubChemNumberPatents"],
               "LocalDatabasePath": "/path/to/db",
               "ResultsPath": "/path/to/results",
               "NumberThreads": 4,
               "PeakListPath": "/path/to/peaks"}
        )

        self.assertIn("FIO00665_peaks.csv", out)
        self.assertIn("FIO00665_config.txt", out)
        self.assertIn("PeakListPath=/path/to/peaks/FIO00665_peaks.csv", out["FIO00665_config.txt"])
        self.assertIn("MetFragScoreWeights=0.8,0.2\n", out["FIO00665_config.txt"])
        self.assertIn("MetFragScoreTypes=FragmenterScore,PubChemNumberPatents\n", out["FIO00665_config.txt"])
        self.assertIn("PrecursorIonMode=-1\n", out["FIO00665_config.txt"])
        self.assertIn("IsPositiveIonMode=False\n", out["FIO00665_config.txt"])

        # Spectrum 2 --------------------
        out = MBSpectrum("./example_massbank_records/EQ308406.txt")._to_metfrag_format(
            **{"MetFragScoreWeights": [1.0],
               "MetFragScoreTypes": ["FragmenterScore"],
               "LocalDatabasePath": "/path/to/db",
               "ResultsPath": "/path/to/results",
               "NumberThreads": 4,
               "PeakListPath": "/path/to/peaks"}
        )

        self.assertIn("EQ308406_peaks.csv", out)
        self.assertIn("EQ308406_config.txt", out)
        self.assertIn("PeakListPath=/path/to/peaks/EQ308406_peaks.csv", out["EQ308406_config.txt"])
        self.assertIn("MetFragScoreWeights=1.0\n", out["EQ308406_config.txt"])
        self.assertIn("MetFragScoreTypes=FragmenterScore\n", out["EQ308406_config.txt"])
        self.assertIn("PrecursorIonMode=1\n", out["EQ308406_config.txt"])
        self.assertIn("IsPositiveIonMode=True\n", out["EQ308406_config.txt"])

        # Spectrum 3 --------------------
        spectra = []
        acc = []
        for mb_fn in glob.iglob(os.path.join("example_massbank_records", "EA0004[01][0-9].txt")):
            spectra.append(MBSpectrum(mb_fn))
            acc.append(spectra[-1].get("accession"))

        merged_spectrum = MBSpectrum.merge_spectra(spectra)
        out = merged_spectrum._to_metfrag_format(
            **{"MetFragScoreWeights": [1.0],
               "MetFragScoreTypes": ["FragmenterScore"],
               "LocalDatabasePath": "/path/to/db",
               "ResultsPath": "/path/to/results",
               "NumberThreads": 4,
               "PeakListPath": "/path/to/peaks"}
        )

        peaks_fn = merged_spectrum.get("accession") + "_peaks.csv"
        config_fn = merged_spectrum.get("accession") + "_config.txt"

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
        rt_ref = [spec.get("retention_time")]
        precmz_ref = [spec.get("precursor_mz")]
        recordtitle_ref = [spec.get("record_title")]
        ce_ref = [spec.get("collision_energy")]

        # -------------------
        # WITH RT AGGREGATION
        # -------------------
        merged_spectrum = MBSpectrum.merge_spectra(spectra, rt_agg_fun=np.min)  # type: MBSpectrum

        self.assertEqual("FBZONXHGGPHHIY-UHFFFAOYSA-N", merged_spectrum.get("inchikey"))
        self.assertEqual(acc_ref, merged_spectrum.get("original_accessions"))
        self.assertEqual("FIO", merged_spectrum.get("accession")[:3])
        self.assertEqual(MBSpectrum._get_new_accession_id(merged_spectrum.get("original_accessions")),
                         merged_spectrum.get("accession"))
        self.assertEqual(None, merged_spectrum.get("retention_time_unit"))
        self.assertEqual(None, merged_spectrum.get("retention_time"))
        self.assertEqual(precmz_ref[0], merged_spectrum.get("precursor_mz"))
        self.assertEqual(recordtitle_ref[0], merged_spectrum.get("record_title"))
        self.assertEqual(ce_ref[0], merged_spectrum.get("collision_energy"))

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
        merged_spectrum = MBSpectrum.merge_spectra(spectra, rt_agg_fun=np.min)  # type: MBSpectrum

        self.assertEqual("OUSYWCQYMPDAEO-UHFFFAOYSA-N", merged_spectrum.get("inchikey"))
        self.assertEqual(acc_ref, merged_spectrum.get("original_accessions"))
        self.assertEqual("EA", merged_spectrum.get("accession")[:2])
        self.assertEqual(MBSpectrum._get_new_accession_id(merged_spectrum.get("original_accessions")),
                         merged_spectrum.get("accession"))
        self.assertEqual("min", merged_spectrum.get("retention_time_unit"))
        self.assertEqual(np.min(rt_ref), merged_spectrum.get("retention_time"))
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


if __name__ == '__main__':
    unittest.main()
