import json
import re
import hashlib
import pandas as pd


def _entry_has_rt(entry):
    for meta_info in entry["metaData"]:
        if meta_info["name"] == "retention time":
            return True

    return False


def _entry_is_LC(entry):
    for meta_info in entry["metaData"]:
        if meta_info["name"] == "instrument type" and re.search("LC", meta_info["value"]):
            return True

    return False


if __name__ == "__main__":
    fn = "/run/media/bach/EVO500GB/data/MassBank_from_MONA/MoNA-export-MassBank.json"

    with open(fn, "r") as monafile:
        monadb = json.load(monafile)

    rows = []
    for idx, entry in enumerate(monadb):
        print("%d/%d" % (idx + 1, len(monadb)))

        # Check whether the entry has a retention time
        if not _entry_has_rt(entry):
            continue

        # Check whether the system is LC
        if not _entry_is_LC(entry):
            continue

        # Get entry ID string
        idstr = re.search("[A-Z]{2,3}", entry["id"]).group(0)

        # Get InChIKey
        inchikey = entry["compound"][0]["inchiKey"]

        # Column
        column = None
        for meta_info in entry["metaData"]:
            if meta_info["name"] == "column":
                column = meta_info["value"]
        if not column:
            continue

        column_hash = hashlib.md5(column.replace(" ", "").lower().encode()).hexdigest()[:8]

        rows += [[idstr, inchikey, column_hash, column, entry["id"]]]

    df = pd.DataFrame(rows, columns=["acquisition", "inchikey", "column_hash", "column", "entryid"])

    print(df.head())

    df.to_csv("stats.csv", index=False)
