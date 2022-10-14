from rdkit import Chem
from rdkit.Chem import rdfiltercatalog

def load_functional_groups():
    """
    Loads the functional groups from a file.
    """
    QueryGroups = rdfiltercatalog.GetFlattenedFunctionalGroupHierarchy(normalized=False)
    flattened_functional_groups = list(QueryGroups.items())

    new_functional_groups = [
                            ("Halogen.Iodine.Aromatic", Chem.MolFromSmarts('[I;$(I-!@c)]')),
                            ("Halogen.Fluorine.Aromatic", Chem.MolFromSmarts('[F;$(F-!@c)]')),
                            ("Halogen.Chlorine.Aromatic", Chem.MolFromSmarts('[Cl;$(Cl-!@c)]')), 
                            ("Nitrile", Chem.MolFromSmarts('[$(C-!@[#6])](#N)')),
                            ("Nitrile.Aromatic", Chem.MolFromSmarts('[$(C-!@c)](#N)')),
                            ("Nitrile.Aliphatic", Chem.MolFromSmarts('[$(C-!@C)](#N)')),
                            ('Acylcyanide', Chem.MolFromSmarts('[$([#6]-!@[#6])](=O)C#N')),
                            ('Acylcyanide.Aromatic', Chem.MolFromSmarts('[$([#6]-!@c)](=O)C#N')),
                            ('Acylcyanide.Aliphatic', Chem.MolFromSmarts('[$([#6]-!@C)](=O)C#N')),
                            ('Ketone', Chem.MolFromSmarts('[$(C-!@[#6])](=O)[#6]')),
                            ('Ketone.Aromatic', Chem.MolFromSmarts('[$(C-!@[a])](=O)[#6]')),
                            ('Ketone.Aliphatic', Chem.MolFromSmarts('[$(C-!@[C])](=O)[#6]')),
                            ('Thiol', Chem.MolFromSmarts('[S;H1;$(O-!@[#6;!$(C=!@[O,N])])]')),
                            ('Thiol.Aromatic', Chem.MolFromSmarts('[S;H1;$(S-!@c)]') ),
                            ('Thiol.Aliphatic', Chem.MolFromSmarts('[S;H1;$(O-!@[C;!$(C=!@[O,N])])]')),
                            ('Thioester', Chem.MolFromSmarts('S([#6])[$(C-!@[#6])](=O)')),
                            ('Thioester.Aromatic', Chem.MolFromSmarts('S([#6])[$(C-!@c)](=O)')),
                            ('Thioester.Aliphatic', Chem.MolFromSmarts('S([#6])[$(C-!@C)](=O)')),
                            ('Ether', Chem.MolFromSmarts('[$(O-!@[#6])]([#6])[#6]')),
                            ('Ether.Aromatic', Chem.MolFromSmarts('[$(O-!@c)]([#6])[#6]')),
                            ('Ether.Aliphatic', Chem.MolFromSmarts('[$(O-!@C)]([#6])[#6]')),
                            ('Ester', Chem.MolFromSmarts('O([#6])[$(C-!@[#6])](=O)')),
                            ('Ester.Aromatic', Chem.MolFromSmarts('O([#6])[$(C-!@c)](=O)')),
                            ('Ester.Aliphatic', Chem.MolFromSmarts('O([#6])[$(C-!@[C])](=O)')),
                            ('Phenyl', Chem.MolFromSmarts('c1cccc[$(c-!@[#6])]1')),
                            ('Phenyl.Aromatic', Chem.MolFromSmarts('c1cccc[$(c-!@c)]1')),
                            ('Phenyl.Aliphatic', Chem.MolFromSmarts('c1cccc[$(c-!@C)]1')),
                            ('Amide', Chem.MolFromSmarts('[OX1]=[$(C-!@[#6])][$(N-!@[#6])]')),
                            ('Amide.Aromatic', Chem.MolFromSmarts('[OX1]=[$(C-!@c)][$(N-!@c)]')),
                            ('Amide.Aliphatic', Chem.MolFromSmarts('[OX1]=[$(C-!@C)][$(N-!@C)]')),
                            ('Alkene', Chem.MolFromSmarts('[$(C-!@[#6])](=C)')),
                            ('Alkene.Aromatic', Chem.MolFromSmarts('[$(C-!@c)](=C)')),
                            ('Alkene.Aliphatic', Chem.MolFromSmarts('[$(C-!@C)](=C)')),
                            ('Alkyne', Chem.MolFromSmarts('[$(C-!@[#6])](#C)')),
                            ('Alkyne.Aromatic', Chem.MolFromSmarts('[$(C-!@c)](#C)')),
                            ('Alkyne.Aliphatic', Chem.MolFromSmarts('[$(C-!@C)](#C)')),
                            ('Phosphonic acid', Chem.MolFromSmarts('[$(P-!@[#6])]([#8])([#8])=O')),
                            ('Phosphonic acid.Aromatic', Chem.MolFromSmarts('[$(P-!@c)]([#8])([#8])=O')),
                            ('Phosphonic acid.Aliphatic', Chem.MolFromSmarts('[$(P-!@C)]([#8])([#8])=O')),
                            ('Phosphate', Chem.MolFromSmarts('[$(O-!@[#6])]P([#8])([#8])=O')),
                            ('Phosphate.Aromatic', Chem.MolFromSmarts('[$(O-!@c)]P([#8])([#8])=O')),
                            ('Phosphate.Aliphatic', Chem.MolFromSmarts('[$(O-!@C)]P([#8])([#8])=O')),
                            ('Phosphodiester', Chem.MolFromSmarts('[$(O-!@[#6])]P([#8][#6])([#8])=O')),
                            ('Phosphodiester.Aromatic', Chem.MolFromSmarts('[$(O-!@c)]P([#8][#6])([#8])=O')),
                            ('Phosphodiester.Aliphatic', Chem.MolFromSmarts('[$(O-!@C)]P([#8][#6])([#8])=O')),
                            ('Sulfoxide', Chem.MolFromSmarts('[$([#16X3]-!@[#6])]=[OX1]')),
                            ('Sulfoxide.Aromatic', Chem.MolFromSmarts('[$([#16X3]-!@c)]=[OX1]')),
                            ('Sulfoxide.Aliphatic', Chem.MolFromSmarts('[$([#16X3]-!@C)]=[OX1]')),
                            ('Sulfone', Chem.MolFromSmarts('[$([#16X4]-!@[#6])](=[OX1])=[OX1]')), #Sulfone. Low specificity: Hits all sulfones, including heteroatom-substituted sulfones: sulfonic acid, sulfonate, sulfuric acid mono- & di- esters, sulfamic acid, sulfamate, sulfonamide
                            ('Sulfone.Aromatic', Chem.MolFromSmarts('[$([#16X4]-!@c)](=[OX1])=[OX1]')),
                            ('Sulfone.Aliphatic', Chem.MolFromSmarts('[$([#16X4]-!@C)](=[OX1])=[OX1]')),
                            ('Haloalkyl', Chem.MolFromSmarts('F[$(C-!@[#6])](F)(F)')),
                            ('Haloalkyl.Aromatic', Chem.MolFromSmarts('F[$(C-!@c)](F)(F)')),
                            ('Haloalkyl.Aliphatic', Chem.MolFromSmarts('F[$(C-!@C)](F)(F)'))                             
    ]

    flattened_functional_groups.extend(new_functional_groups)
    flattened_functional_groups = sorted(flattened_functional_groups)

    return flattened_functional_groups


def get_single_omp_FGv3(smiles, func_smarts='[F,Cl,Br,I]', group_catalog=load_functional_groups(), display_mol=False):
    """
    #https://github.com/rdkit/rdkit/blob/4fafb72212163bbdaeafbb364d7f922b95a823bf/Code/GraphMol/FilterCatalog/Wrap/rough_test.py

    Returns a dictionary of functional groups and matches for ortho/meta/para
    based on the funct_smarts.
    Example of usage:
    get_single_omp_FG(regioselective_dataframe['Halide'][0], display_mol=False) -->
    {'Halogen.Aromatic': [('ortho', (6, 5, 7, 8)), ('meta', (8, 7, 9, 10, 11)), ('para', (6, 5, 4, 3, 10, 11))], 
    'Halogen.NotFluorine.Aromatic': [('ortho', (6, 5, 7, 8)), ('meta', (8, 7, 9, 10, 11)), ('para', (6, 5, 4, 3, 10, 11))],
    'Nitro.Aromatic': [('ortho', (11, 10, 3, 1, 0, 2)), ('meta', (6, 5, 4, 3, 1, 0, 2)), ('para', (8, 7, 5, 4, 3, 1, 0, 2))]}
    """

    mol = Chem.MolFromSmiles(smiles)

    hits_FG = [(name, Chem.MolToSmarts(pat)) for name, pat in group_catalog if mol.HasSubstructMatch(pat)]
    #print('hits_FG')
    #print(hits_FG)
    
    # hits_FG_getsubstructurematches = [(name, Chem.MolToSmarts(pat)) for name, pat in group_catalog if mol.GetSubstructMatches(pat)]
    # print('hits_FG_getsubstructurematches')
    # print(hits_FG_getsubstructurematches)
    FG_group_o_m_p = dict()

    for group, smarts in hits_FG:
        o_m_p_list = []
        #print(group, smarts)
        ortho_smarts = Chem.MolFromSmarts(func_smarts+'cc'+smarts)
        meta_smarts = Chem.MolFromSmarts(func_smarts+'ccc'+smarts)
        para_smarts = Chem.MolFromSmarts(func_smarts+'cccc'+smarts)
        #print(f"ortho smart: {Chem.MolToSmarts(ortho_smarts)}")
        
        ortho_match = mol.GetSubstructMatches(ortho_smarts)
        if len(ortho_match) != 0:
            o_m_p_list.append(('ortho', ortho_match))
            #print(f'{group} ortho match at: {ortho_match}')
            if display_mol:
                display(mol)
        meta_match = mol.GetSubstructMatches(meta_smarts)
        if len(meta_match) != 0:
            o_m_p_list.append(('meta', meta_match))
            #print(f'{group} meta match at: {meta_match}')
            if display_mol:
                display(mol)
        #para only needs one match else it will count the same match twice - mol.GetSubstructMatches(para_smarts)        
        para_match = mol.GetSubstructMatch(para_smarts)
        if len(para_match) != 0:
            o_m_p_list.append(('para', para_match))
            #print(f'{group} para match at: {para_match}')
            if display_mol:
            #display(SVG(draw_chemical_reaction(col[0],  highlightByReactant=True)))
                display(mol)
        FG_group_o_m_p[group] = o_m_p_list

    return FG_group_o_m_p