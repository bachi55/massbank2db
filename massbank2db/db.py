import sqlite3

def create_db(file_pth):
    """

    :return:
    """
    conn = sqlite3.connect(file_pth)

    # Molecules Table
    conn.execute(
        "CREATE TABLE molecules( \
             cid               INTEGER PRIMARY KEY, \
             inchi             VARCHAR NOT NULL, \
             inchikey          VARCHAR NOT NULL, \
             inchikey1         VARCHAR NOT NULL, \
             inchikey2         VARCHAR NOT NULL, \
             inchikey3         VARCHAR NOT NULL, \
             smiles_iso        VARCHAR NOT NULL, \
             smiles_can        VARCHAR NOT NULL, \
             molecular_weight  FLOAT NOT NULL, \
             exact_mass        FLOAT NOT NULL, \
             molecular_formula VARCHAR NOT NULL, \
             )")
    conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey_index ON molecules(inchikey)")
    conn.execute("CREATE INDEX IF NOT EXISTS molecules_inchikey1_index ON molecules(inchikey1)")
    conn.execute("CREATE INDEX IF NOT EXISTS molecules_exact_mass_index ON molecules(exact_mass)")
    conn.execute("CREATE INDEX IF NOT EXISTS molecules_mf_index ON molecules(molecular_formula)")

    # Datasets Meta-information table
    conn.execute(
        "CREATE TABLE datasets( \
            name                VARCHAR PRIMARY KEY, \
            contributor         VARCHAR NOT NULL, \
            retention_time_unit VARCHAR NOT NULL, \
            copyright           VARCHAR, \
            license             VARCHAR, \
            column_name         VARCHAR, \
            column_temperature  VARCHAR, \
            flow_gradient       VARCHAR, \
            flow_rate           VARCHAR, \
            solvent_A           VARCHAR, \
            solvent_B           VARCHAR, \
            solvent             VARCHAR, \
            instrument_type     VARCHAR, \
            instrument          VARCHAR)"
    )

    # Spectra Meta-information table
    conn.execute(
        "CREATE TABLE spectra_meta( \
            accession           VARCHAR PRIMARY KEY, \
            dataset             VARCHAR NOT NULL, \
            record_title        VARCHAR NOT NULL, \
            molecule            INTEGER NOT NULL, \
            ion_mode            VARCHAR NOT NULL, \
            precursor_mz        FLOAT NOT NULL, \
            precursor_type      FLOAT NOT NULL, \
            retention_time      FLOAT NOT NULL, \
            collision_energy    FLOAT, \
            ms_level            INTEGER, \
            resolution          FLOAT, \
            fragmentation_type  VARCHAR, \
            origin              VARCHAR, \
         FOREIGN KEY(molecule) REFERENCES molecules(cid),\
         FOREIGN KEY(dataset) REFERENCES datasets(name))"
    )
