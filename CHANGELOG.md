# Changelog

## Version 0.6.2
**Minor changes**

### Spectra merging
- intensities can be normalized (to maximum one) *after* the spectra have been merged
- 

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
