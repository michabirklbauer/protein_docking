#!/usr/bin/env python3

# CREATE DATASET FROM PDB ENTRIES
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

def create_dataframe(input_name, output_name, mol_name, EC_nr):

    with open(input_name, "r") as f:
        content = f.read()
        f.close()

    fields = content.split(";")
    if fields[5].split(":")[0].strip() == "PDB_ENTRIES":
        entries = fields[5].split(":")[1].split(",")
    else:
        print("Wrong field selected for entries!")
        print(fields[5].split(":")[0])
        return 1

    if fields[2].split(":")[0].strip() == "PDB_ORGANISM":
        organism = fields[2].split(":")[1].strip()
    else:
        print("Wrong field selected for organism!")
        print(fields[2].split(":")[0])
        return 1

    INDEX = []
    NAME = []
    EC_NR = []
    PDB_ENTRIES = []
    PDB_ORGANISM = []

    i = 0
    for entry in entries:
        INDEX.append(i)
        i = i + 1
        NAME.append(mol_name)
        EC_NR.append(EC_nr)
        PDB_ENTRIES.append(entry.strip())
        PDB_ORGANISM.append(organism)

    df = pd.DataFrame(
            {"INDEX": INDEX,
             "NAME": NAME,
             "EC_NR": EC_NR,
             "PDB_ENTRY": PDB_ENTRIES,
             "PDB_ORGANISM": PDB_ORGANISM}
        )

    train, test = train_test_split(df, train_size = 0.8, random_state = 42, shuffle = True)

    train.to_csv(output_name + "_train.csv", index = False)
    test.to_csv(output_name + "_test.csv", index = False)

    return 0

def create_dataframe_human(input_name, output_name, mol_name, EC_nr):

    with open(input_name, "r") as f:
        content = f.read()
        f.close()

    fields = content.split(";")
    if fields[5].split(":")[0].strip() == "PDB_ENTRIES":
        entries = fields[5].split(":")[1].split(",")
    else:
        print("Wrong field selected for entries!")
        print(fields[5].split(":")[0])
        return 1

    if fields[2].split(":")[0].strip() == "PDB_ORGANISM":
        organism = fields[2].split(":")[1].strip()
    else:
        print("Wrong field selected for organism!")
        print(fields[2].split(":")[0])
        return 1

    INDEX = []
    NAME = []
    EC_NR = []
    PDB_ENTRIES = []
    PDB_ORGANISM = []

    i = 0
    for entry in entries:
        INDEX.append(i)
        i = i + 1
        NAME.append(mol_name)
        EC_NR.append(EC_nr)
        PDB_ENTRIES.append(entry.strip())
        PDB_ORGANISM.append(organism)

    df = pd.DataFrame(
            {"INDEX": INDEX,
             "NAME": NAME,
             "EC_NR": EC_NR,
             "PDB_ENTRY": PDB_ENTRIES,
             "PDB_ORGANISM": PDB_ORGANISM}
        )

    df.to_csv(output_name + "_human.csv", index = False)

    return 0

if __name__ == "__main__":
    r = create_dataframe("PDB2.txt", "data", "Cyclooxygenase-2", "1.14.99.1")
    print(r)
    r = create_dataframe_human("PDB.txt", "data", "Cyclooxygenase-2", "1.14.99.1")
    print(r)
