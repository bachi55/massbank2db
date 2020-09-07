# Database Layout

## Dataset Characterizations

The different [contributors](https://github.com/MassBank/MassBank-data/blob/main/List_of_Contributors_Prefixes_and_Projects.md) and their prefix ids are separated into different datasets. Furthermore, the entries of a particular contributor are grouped by their measurement setup:
```sqlite
SELECT COUNT(DISTINCT inchikey), GROUP_CONCAT(accession) FROM information
   GROUP BY ion_mode, instrument, instrument_type, column_name, column_temperature, 
            flow_gradient, flow_rate, solvent_A, solvent_B
```
Aim of this separation is, to ensure that the (tandem) mass spectra (MS) *and* retention times (RT) are measured under the same conditions. For the RTs this is particularly important, as RTs for the same compounds can be totally different in different chromatographic setups. Please note, the separation into different datasets does not impose any restrictions, but rather makes it possible to query (MS, RT)-tuples measured under the same conditions, if this is somehow useful in the application.

## Tables

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

### Table: datasets

As previously described, the Massbank entries are grouped by contributors and measurement setups. The tabel *datasets* stores meta-information to all those groups.

| Column | Description | Primary Key | Foreign Key | Index |
