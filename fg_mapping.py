#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to detect functional groups in molecules using RDKit functionality and defined functional group SMARTS in a
separate file (`Functional_Group_Hierarchy.txt`).

Authors: Alex Müller
Last Update: 2022-08-09
"""
import os
import argparse

import pandas as pd
from rdkit.Chem import MolFromSmiles, ParseMolQueryDefFile, CanonSmiles


def read_fg_hierarchy(fg_file: str = "./Functional_Group_Hierarchy.txt"):
    """
    Load and process the functional group hierarchy definition

    :param fg_file: Filename of the file containing the functional group hierarchy definition.
    :return: {pandas.DataFrame} DataFrame containing names and patterns of the functional groups
    """
    # read the functional group hierarch file using RDKit's parser and pandas
    fg_dict = ParseMolQueryDefFile(fg_file)
    fg = pd.DataFrame({"Name": fg_dict.keys(), "Pattern": fg_dict.values()})
    fg["Level"] = fg["Name"].str.count(r'\.')
    return fg.sort_values("Name")


def match_fg_hierarchy(smiles: str, hierarchy: pd.DataFrame, lowest: bool = True, atoms: bool = False):
    """
    Match all functional groups in `hierarchy` on the given molecule

    :param smiles: {str} a SMILES string for one molecule
    :param hierarchy: {pandas.DataFrame} a DataFrame containing the functional group hierarchy
    :param lowest: {bool} only return the lowest level matches
    :param atoms: {bool} whether atom indices for the match should be returned
    :return: {dict} dictionary containing names of the matches and corresponding counts
    """
    low_match = {}  # dict to store only the lowest hierarchy match per atom
    all_match = {}  # dict to count and store all matches
    try:
        m = MolFromSmiles(CanonSmiles(smiles))
    except TypeError:
        return {None: 0}
    for i, row in hierarchy.iterrows():
        match = m.GetSubstructMatches(row["Pattern"])
        if len(match):
            for mtch in match:
                for atm_idx in mtch:
                    low_match[atm_idx] = row["Name"]  # sort tuple of atom matches to match lower levels
                all_match[row["Name"]] = len(match)  # store count
    if lowest:
        if atoms:
            return {v: (all_match[v], sorted([k2 for k2, v2 in low_match.items() if v2 == v])) for v in
                    low_match.values()}
        return {v: all_match[v] for v in low_match.values()}
    return all_match


def smiles_fg_counts(smiles: list, fg_file: str = "./Functional_Group_Hierarchy.txt", lowest: bool = True,
                     atoms: bool = False):
    """
    Matches functional groups in a list of molecules represented as SMILES strings.

    :param smiles: {list} list of SMILES strings for which functional groups should be matched
    :param fg_file: Filename of the file containing the functional group hierarchy definition.
    :param lowest: {bool} only return the lowest level matches
    :param atoms: {bool} whether atom indices for the match should be returned
    :return: {dict} Dictionary with functional group matches for every SMILES string
    """
    fg_hrchy = read_fg_hierarchy(fg_file)
    return {str(s): match_fg_hierarchy(s, fg_hrchy, lowest, atoms) for s in smiles}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fragment Mapping')
    parser.add_argument('-f', '--fg_def', type=str, default="./Functional_Group_Hierarchy.txt",
                        help='Functional group hierarchy definition filename')
    parser.add_argument('-i', '--input', type=str, required=True, help='Smiles input filename with the headers '
                                                                       '"Index" and "SMILES"')
    parser.add_argument('-s', '--sep', type=str, default='\t', help="Column separator in input file")
    parser.add_argument('--start', type=int, default=1, help="Index of the first entry to process")
    parser.add_argument('--stop', type=int, default=None, help="Index of the last entry to process, None=all")
    parser.add_argument('-o', '--output', type=str, required=True, help='Functional group mapping output filename')
    args = parser.parse_args()

    if not args.fg_def.startswith('/'):
        fg_def = os.path.join(os.path.dirname(__file__), args.fg_def)
    else:
        fg_def = args.fg_def

    data = pd.read_csv(args.input, delimiter=args.sep, index_col=0)
    if args.stop:
        data = data.loc[args.start:args.stop]

    matches = smiles_fg_counts(data["SMILES"].tolist(), fg_def, True, False)

    with open(args.output, "w") as outf:
        outf.write('Index\tSMILES\tFG_Name\tFG_Count\n')
        for idx, smls in data["SMILES"].iteritems():
            for fgrp, cnt in matches[str(smls)].items():
                if fgrp:
                    outf.write(str(idx) + '\t' + str(smls) + '\t' + fgrp + '\t' + str(cnt) + '\n')
