# Build a local SQLite Database from Massbank

[Massbank](https://github.com/MassBank/MassBank-data) is an open repository for mass spectra (MS) including tandem MS (MS/MS). Meta information, such as measured compound, retention times (RT) and instrument setup, are provided with the MS. 

**massbank2db** provides functionality to parse Massbank [entries](https://github.com/MassBank/MassBank-data/blob/main/CASMI_2016/SM800003.txt), filter them (e.g. removing entries without RT) and organize them into a local SQLite database. This project took inspiration from the similar project [msp2db](https://github.com/computational-metabolomics/msp2db). However, **massbank2db** differs in some main-features:

  - DB stores more meta information
  - Data filtering, e.g. only include complete entries with MS2 and RT
  - Load the compound information, such as InChI and SMILES, from ([a local copy](https://github.com/bachi55/local_pubchem_db) of) PubChem using the CID's of the Massbank entries
  - Add molecular candidates from PubChem for each MS or MS/MS spectrum
