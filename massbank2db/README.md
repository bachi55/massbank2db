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
| --- | --- | --- | --- | --- |
| name | Dataset identifier compiled from the contributors' accession-prefixes and a running sub-group integer, e.g. AU_001 | True | None | True | 
| contributor | Contributor | False | None | False | 
| ion_mode | Ionization mode. Either "positive" or "negative" | False | None | False | 
| num_spectra | Number of MS(/MS) spectra in this dataset | False | None | False | 
| num_compounds | Number of unqiue compounds in this dataset | False | None | False | 
| copyright | [Massbank copyright information](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#2.1.6) | False | None | False | 
| license | [Massbank license information](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#215-license) | False | None | False | 
| column_name | [Commercial Name of Chromatography Column and Manufacture](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#246-subtag-column_name) | False | None | False | 
| column_type | Either HILIC or Reversed Phased (RP) (TODO: Information needs to be extracted manualy) | False | None | False | 
| column_temperature | [Static Column Temperature](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#246-subtag-column_temperature) | False | None | False | 
| flow_gradient | [Gradient of Mobile Phases in LC-MS](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#246-subtag-flow_gradient) | False | None | False |
| flow_rate | [Flow Rate of Migration Phase](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#246-subtag-flow_rate) | False | None | False | 
| solvent_A | [Chemical Composition of Buffer Solution (A)](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#246-subtag-solvent) | False | None | False | 
| solvent_B | Chemical Composition of Buffer Solution (B) | False | None | False | 
| solvent | Chemical Composition of Buffer Solution (if no separate information for A and B are given) | False | None | False | 
| instrument_type | [General Type of Instrument](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#242-acinstrument_type) | False | None | 
| instrument | [Commercial Name and Model of Chromatographic Separation Instrument, if any were coupled, and Mass Spectrometer and Manufacturer](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#241-acinstrument) | False | None | False |
| column_dead_time | Estimated column dead-time (TODO: calculate from LC specification or estimate) | False | None | False | 

### Table: spectra_meta

For each MS(/MS) spectrum the DB stores meta-information, such as collision enegery or precursor m/z, as well a reference to the corresponding molecule. 

| Column | Description | Primary Key | Foreign Key | Index |
| --- | --- | --- | --- | --- |
| accession | Massbank ID of the spectrum file | True | None | True | 
| dataset | Dataset identifier | False | datasets(name) | False | 
| record_title | [Brief Description of Massbank entry](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#2.1.2) | False | None | False | 
| molecule | CID of the molecule associated with the entry | False | molecules(cid) | False | 
| precursor_mz | [m/z of Precursor Ion in MSn spectrum](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#251-subtag-precursor_mz) | False | None | False | 
| precursor_type | [Type of Precursor Ion in MSn spectrum](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#251-subtag-precursor_type) | False | None | False | 
| collision_energy | [Collision Energy for Dissociation](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#245-subtag-collision_energy) | False | None | False | 
| ms_type | [MS(1) no fragmentation, MS2 1st generation product ion spectrum, ...](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#243-acmass_spectrometry-ms_type) | False | None | False |
| resolution | [Resolution (aka Mass Resolution or Resolving Power) is the smallest mass difference between two equal magnitude peaks so that the valley between them is a specified fraction of the peak height](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#245-subtag-resolution) | False | None | False |
| fragmentation_type | [Fragmentation method used for dissociation or fragmentation](https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md#245-subtag-fragmentation_mode) | False | None | False | 

### Table: spectra_peaks

This tables stores all peaks, i.e. (m/z, int)-tuples, of Massbank and their associated accessions.

| Column | Description | Primary Key | Foreign Key | Index |
| --- | --- | --- | --- | --- | 
| spectrum | Spectrum id | False | spectra_meta(accession) | True | 
| mz | Mass per charge (m/z) of the peak | False | None | False | 
| itensity | Intensiry of the peak (unnormalized) | False | None | False | 

### Table: spectra_raw_rts

This tables stores the retention times associated with the Massbank entries.

| Column | Description | Primary Key | Foreign Key | Index |
| --- | --- | --- | --- | --- | 
| spectrum | Spectrum id | False | spectra_meta(accession) | True | 
| retention_time | Retention time | False | None | True | 
| retention_time_unit | Retention time unit, e.g. min or sec | False | None | False | 

### Table: spectra_candidates

This table associates each spectrum in the DB with a set of molecular candidate structures extracted from PubChem. The candidate sets ```C_i``` for spectrum ```i``` are comprised of all molecular structurs those exact mass ```m_is``` is within an +/- X ppm window around the (estimated) exact mass ```m_i = prec_i - ion_mode_i * adduct_mass_i``` of the molecule associated with the spectrum: ```m_is € [m_i - (m_i * X) / 1e6, m_i + (m_i * X) / 1e6]``` (following the MetFragRelaunched library, [here](https://github.com/ipb-halle/MetFragRelaunched/blob/c57f9d2b406350b2357ce9f7ce42a286cefcca13/MetFragLib/src/main/java/de/ipbhalle/metfraglib/additionals/MathTools.java#L16) and [here](https://github.com/ipb-halle/MetFragRelaunched/blob/c57f9d2b406350b2357ce9f7ce42a286cefcca13/MetFragLib/src/main/java/de/ipbhalle/metfraglib/database/LocalMySQLDatabase.java#L23)).

| Column | Description | Primary Key | Foreign Key | Index |
| --- | --- | --- | --- | --- | 
| spectrum | Spectrum id | False | spectra_meta(accession) | True | 
| candidate | CID of the molecule associated with the entry | False | molecules(cid) | True | 
| mf_gt_equal | If 1 (=True), the candidate's molecular formula is equal to the *ground truth* molecular formula of the molecule associated with the spectrum | False | None | False | 
| mf_pred_equal | If 1 (=True), the candidate's molecular formula is equal to the *predicted* molecular formula of the molecule associated with the spectrum | False | None | False | 
| ppm_diff_gt | Mass difference (in ppm) of the candidate and the *ground truth* exact mass of the spectrum | False | None | True |
| ppm_diff_est | Mass difference (in ppm) of the candidate and the *estimated* exact mass of the spectrum | False | None | True |
