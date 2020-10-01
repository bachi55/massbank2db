# Changelog

## Version 0.4.2
- fix bug in new accession id generator: Total length now fixed to 8, also when acc. pref. has length 3.
- dataset generation now takes the 'fragmentation_mode' into account
- fix some broken tests

## Version 0.4.1
- modify spectra grouping small 
- small fix to the pubchem information retrieval when iterating the spectra
- changes in the metfrag output format

## Version 0.4.0
- simplify the MassbankDB class constructor parameters
- 'insert_dataset' and 'iter_spectra' now have these parameters

## Version 0.3.0
**Feature release**
- add support to output MetFrag format for spectra and candidates

## Version 0.2.0
**Feature release**
- adressed [#2](https://github.com/bachi55/massbank2db/issues/2): Implement a spectra merging strategy

## Version 0.1.1
- lot of minor improvements
- added tutorials: "How to build a DB?" and "How to access data from the DB?"

## Version 0.1.0
**Initial release**
- implements DB building
- implements spectra query 
