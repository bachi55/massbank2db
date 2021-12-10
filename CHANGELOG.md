# Changelog

## Version 0.8.0
**Feature release**

- Add importer for CFM-ID insilico spectra to MBSpectrum objects

## Version 0.7.1

- Add support for some new adduct types for the MetFrag tool exporter (based on [this source](https://github.com/ipb-halle/MetFragRelaunched/blob/a61900b82c0cd5835f0b7d6ae822054b1e3eb104/MetFragTools/src/main/java/de/ipbhalle/metfrag/conversion/ConvertMGFtoMetFragRecord.java)) 

## Version 0.7.0
**Feature release**

### Add support for the CFM-ID in-silico tool

- spectra can be exported to fit the input format of CFM-ID
- multiple collision energies in Massbank are merged and repeated for the three energy levels that CFM-ID expects
- candidate lists can be created

## Version 0.6.4

- Change License from MIT to GNU GPLv3

## Version 0.6.3
**Bug fix release**

- When writing MetFrag candidate sets, the InChIKey1 and InChIKey2 column is now stored as well.

## Version 0.6.2
**Minor changes** 

- **Spectra merging** intensities can be normalized (to maximum one) *after* the spectra have been merged

## Version 0.6.1
**Bug fix release**
- fixed several bugs in the SIRIUS exporter
- updated the DB usage tutorial to be compatible with the latest version

## Version 0.6.0
**Feature release**
- add [exporting tool for SIRIUS .ms-files](https://github.com/bachi55/massbank2db/issues/4)
- add more tests

## Version 0.5.0
**Feature release**
- add [pipeline design](https://github.com/bachi55/massbank2db/pull/9) to filter spectra and datasets from the database
- add monoisotopic-mass and xlogp3 (from PubChem) to the compound table
- calculate the mass error between theoretical and measured precursor m/z
- add [estimator for the column dead time](https://github.com/bachi55/massbank2db/pull/10)
- [store all retention times in minutes](https://github.com/bachi55/massbank2db/pull/11) (automatic conversion of RTs in seconds)

## Version 0.4.2
- fix bug in new accession id generator: Total length now fixed to 8, also when acc. pref. has length 3.
- dataset generation now takes the 'fragmentation_mode' into account
- fix some broken tests

## Version 0.4.1
- modify spectra grouping 
- a small fix to the PubChem information retrieval when iterating the spectra
- changes in the MetFrag output format

## Version 0.4.0
- simplify the MassbankDB class constructor parameters
- 'insert_dataset' and 'iter_spectra' now have these parameters

## Version 0.3.0
**Feature release**
- add support to output MetFrag format for spectra and candidates

## Version 0.2.0
**Feature release**
- addressed [#2](https://github.com/bachi55/massbank2db/issues/2): Implement a spectra merging strategy

## Version 0.1.1
- lot of minor improvements
- added tutorials: "How to build a DB?" and "How to access data from the DB?"

## Version 0.1.0
**Initial release**
- implements DB building
- implements spectra query 
