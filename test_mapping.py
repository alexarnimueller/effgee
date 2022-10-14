"""Test script to test a couple of FG matches from the hierarchy definition file that used to make problems."""

import pandas as pd
from fg_mapping import smiles_fg_counts


def test_fg_mapping():
    data = pd.read_csv('test_smiles.txt', delimiter='\t')
    matches = smiles_fg_counts(data["SMILES"].tolist(), "./Functional_Group_Hierarchy.txt", True, False)
    rslt = {s: set(m.keys()) for s, m in matches.items()}
    errors = []
    for _, row in data.iterrows():
        if row.FG not in rslt[row.SMILES]:
            errors.append(f"{row.FG} for {row.SMILES} was not found in {rslt[row.SMILES]}!")
    try:
        assert len(errors) == 0
    except AssertionError:
        for e in errors:
            print(e)
