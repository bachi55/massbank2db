import sqlite3


from massbank2db.utils import named_row_factory


if __name__ == "__main__":
    mb_dbfn = "test_DB.sqlite"

    con = sqlite3.connect(mb_dbfn)
    con.row_factory = named_row_factory

    res = con.execute("SELECT * FROM spectra_meta LIMIT 5")

    for row in res:
        print(row)

    con.row_factory = None
    res = con.execute("SELECT * FROM spectra_meta LIMIT 5")

    for row in res:
        print(row)
