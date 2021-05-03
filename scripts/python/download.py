#!/usr/bin/env python3

# PDB DOWNLOAD SCRIPT
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import urllib.request as ur

def download_pdb_files(filename, output_path = "structures"):

    pdb_link = "https://files.rcsb.org/download/"

    with open(filename, "r") as f:
        content = f.read()
        f.close()

    fields = content.split(";")
    tmp = fields[5]
    if tmp.split(":")[0].strip() == "PDB_ENTRIES":
        PDB_ENTRIES = tmp.split(":")[1].split(",")
    else:
        print("Wrong field selected!")
        print(tmp.split(":")[0])
        return 1

    for entry in PDB_ENTRIES:
        output_name = output_path + "/" + entry.strip() + ".pdb"
        link = pdb_link + entry.strip() + ".pdb"
        ur.urlretrieve(link, output_name)

    return 0

# if __name__ == "__main__":
#     r = download_pdb_files(...)
#     print(r)
