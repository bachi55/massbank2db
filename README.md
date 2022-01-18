# Build a local SQLite Database from Massbank

[MassBank](https://github.com/MassBank/MassBank-data) is an open repository for mass spectrum (MS) data, including 
tandem MS (MS²). Meta information, such as the ground truth molecular structure, retention times (RT) and instrument 
setup, are provided for the MassBank records.   

## The 'massbank2db' package in a nutshell

The **massbank2db** package provides functionality to parse Massbank [records](https://github.com/MassBank/MassBank-data/blob/main/CASMI_2016/SM800003.txt)
and organize them into a local SQLite database. The main focus is to make MassBank's MS data available for machine 
learning method development, which is achieved by the following functionality:

- Molecular structure representations, such as SMILES and InChI, can be updated using PubChem, but matching the 
  provided InChIKeys in MassBank, which allows for consistent and standardized molecule representations for all data 
  points.
- The MassBank records are grouped by their measurement conditions, that is the chromatographic and mass 
  spectrometry setup, such that, for example, retention times are comparable within each group. 
- The database creation process enforces complete information on MassBank records and furthermore filters data which 
  could be corrupted. That means, if relevant information, such as the retention time or ground truth structure 
  annotation is missing, a record is not considered. Furthermore, officially deprecated and otherwise "suspicious" 
  records are filtered out.

## Installation

That's how you install the package:

1) Clone the latest version of the package and change to the directory:
```bash
git clone https://github.com/bachi55/massbank2db
cd massbank2db
```

2) (Optional) Create a **conda** environment:
```bash
conda env create -f environment.yml
conda activate massbank2db
```

3) Install the package:
```bash
# Build the binaries for the spectra merging implementation
python setup.py build_ext --inplace 
# Install the package
pip install .
```

4) (Optional) Test the installation:
```bash
python -m unittest discover -s massbank2db/tests -p 'unittests*.py'

## Expected output ##
# ...s......................
# ----------------------------------------------------------------------
# Ran 26 tests in 0.150s
#
# OK (skipped=1)
```

## Using the package

### Current limitations (<= version 0.9.0)

- The package has been tested with **MassBank release 2020.11** and later releases might cause errors in the record 
  grouping
  - The main reason is that the measurement conditions (specifically for the chromatography) are not reported in a 
    structured and unified style.
  - For the 2020.11 release, a [pull-request](https://github.com/MassBank/MassBank-data/pull/141) homogenising the setup description has been merged.
- The package **does not support MS¹ (only)** records.
  - The main reason is that the package initially was developed for a specific machine learning question ad hand.
  - MS¹ was not relevant for the research question.
  - Adding MS¹ support should not be difficult, but requires some additional thinking.
- The package **requires a local PubChem SQLite DB** with a specific structure. 
  - Again, that is legacy from the original research context. 
  - Such a database can be constructed using the [pubchem2sqlite](https://github.com/bachi55/local_pubchem_db) 
    package using the default DB layout.
- Some **meta-information** for the dataset (= MassBank groups) **cannot be automatically determined**.
  - The liquid-chromatography column type, i.e. reversed-phase (RP) or HILIC, must be added manually. 
  - The column dead-time can only be estimated when all relevant column information are provided. Otherwise, it must 
    be added manually. 
- The package currently **only supports liquid-chromatography (LC)** records.
  - Again, that is legacy from the original research context.
  - Support for gas-chromatography could be added, but it probably requires some tables to be redesigned, i.e. new 
    columns are added. 
  - Specifically the meta-information for the datasets is currently tailored to LC systems.

### Building the MassBank database (example using release 2020.11)

**Note:** At the moment you need to download [forked MassBank release 2020.11](https://github.com/bachi55/MassBank-data/tree/2020.11-branch). 
That is because the fork implements a [small fix in the contributor table](https://github.com/bachi55/MassBank-data/commit/ed2ac2948b97efb1721b05eb89813ffb67107315),
that is allows **masbank2db** to determine all relevant sub-directories in the repository, containing the MS data, 
and associating each sub-directory with the available MassBank record prefixes.

**Note 2:** The following instructions assume, that you have created [your own PubChem SQLite DB](https://github.com/bachi55/local_pubchem_db). 

1) Clone the MassBank repository release 2020.11:
```bash
git clone -b 2020.11-branch https://github.com/bachi55/MassBank-data
```

2) (Optional) Activate your conda environment:
```bash
conda activate massbank2db
```

3) (Optional) Run the data generation script with ```--help``` option:
```bash
python generate_massbank_sqlite.py --help

## Expected output ##  
# usage: generate_massbank_sqlite.py [-h] [--use_pubchem_structure_info] [--pubchem_db_fn PUBCHEM_DB_FN] massbank_repo_dir massbank_db_fn
# 
# positional arguments:
#   massbank_repo_dir     Directory containing the local MassBank copy. It is assumed to be a clone of the MassBank GitHub repository.
#   massbank_db_fn        Output filename of the MassBank SQLite database.
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --use_pubchem_structure_info
#                         Indicating whether a local PubChem copy should be used to update the molecular structure information.
#   --pubchem_db_fn PUBCHEM_DB_FN
#                         Filename of the local PubChem SQLite database file.
```

4) Build the MassBank DB:
```bash
python generate_massbank_sqlite.py /path/to/MassBank-data /tmp/massbank.sqlite --use_pubchem_structure_info --pubchem_db_fn=/path/to/pubchem.sqlite 
```

## Citing the package

If you use this package, please cite:

```bibtex
@software{massbank2db,
  author = {Bach, Eric},
  month = {1},
  title = {{massbank2db: Build a machine learning ready SQLite database from MassBank.}},
  url = {https://github.com/bachi55/massbank2db},
  version = {0.9.0},
  year = {2022}
} 
```