# Build a local SQLite Database from Massbank

[Massbank](https://github.com/MassBank/MassBank-data) is an open repository for mass spectra (MS) including tandem MS (MS/MS). Meta information, such as measured compound, retention times (RT) and instrument setup, are provided with the MS. 

**massbank2db** provides functionality to parse Massbank [entries](https://github.com/MassBank/MassBank-data/blob/main/CASMI_2016/SM800003.txt), filter them (e.g. removing entries without RT) and organize them into a local SQLite database. This project took inspiration from the similar project [msp2db](https://github.com/computational-metabolomics/msp2db). However, **massbank2db** differs in some main-features:
-
