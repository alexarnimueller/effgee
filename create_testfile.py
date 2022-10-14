import os
from rdkit.Chem import SDMolSupplier, MolToSmiles

basedir = r'../EnamineBBs'
f_out = './test_smiles_enamine.txt'

with open(f_out, 'w') as out:
    out.write("SMILES\tFunctionalGroup\n")
    for f in os.listdir(basedir):
        n = f.replace('Enamine_', '').replace('TOP50_', '').split('_')[0]
        print(n)
        mols = SDMolSupplier(os.path.join(basedir, f))

        for i, mol in enumerate([m for m in mols]):
            if i < 50:
                out.write(f"{MolToSmiles(mol, canonical=True)}\t{n}\n")
            else:
                break
