# Database Layout

### Table: molecules 

Each spectrum in Massbank is associated with a single molecule (or compound, terms are used as synonyms). The *molecules* tables stores and entrie for each unique molecules in the Massbank dataset. Thereby, if two spectra are associated with the same molecules only one entry is present, ensuring that the same molecule structures etc. are used. The molecular candidates are stored in this table as well. 

If requested (TODO: add link to code here), all information in this table are extracted from PubChem using the [```CH$LINK: PUBCHEM```](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#2.2.8). Otherwise, the Massbank entry's information is directly used.  

| Column | Description | Primary Key | Foreign Key | Index |
| --- | --- | --- | --- | --- | 
| cid | PubChem compound ID ([PUBCHEM_COMPOUND_CID](https://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_sdtags.pdf)) | True | None | True |
| inchi | InChI representation of the molecule | False | None | False |
| inchikey | Hashed version of the InChI, e.g. UFFBMTHBGFGIHF-UHFFFAOYSA-N  | False | None | True | 
| inchikey1 | First part of the InChIKey, e.g. UFFBMTHBGFGIHF | False | None | True | 
| inchikey2 | Second part of the InChIKey, , e.g. UHFFFAOYSA | False | None | True | 
| smiles_iso | Isomeric SMILES representation | False | None | False |
| smiles_can | Canonical SMILES representation | False | None | False |
| exact_mass | Exact mass of the compound | False | None | True |
| molecular_formula | Molecular formula of the compound | False | None | True |
