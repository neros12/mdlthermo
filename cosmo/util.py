import numpy as np
import os
import pandas as pd
from itertools import chain
from rdkit import Chem
from scipy.linalg import fractional_matrix_power

here = os.path.dirname(os.path.realpath(__file__))

def get_comp_list():
    return pd.read_excel(f'{here}\\model\\component_list.xlsx',
                         sheet_name='component_list', index_col=0)

def get_file(index):
    return f'{here}\\cosmo_file\\{index}.cosmo'

def mgcn_input(smiles, maxatom=60):
    mol = Chem.MolFromSmiles(smiles)
    n_atom = mol.GetNumAtoms()
    assert maxatom >= n_atom, f'The molule has {n_atom} atoms, but the \
maximum number of atoms are {maxatom}.'

    # Node feature matrix
    # atom(11) RS(2) charge(1) bond(4) EZ(2)
    nfm = np.zeros((maxatom, 20))

    nf_dict = {'H':0, 'C':1, 'N':2, 'O':3, 'F':4, # node feature dictionary
               'P':5, 'S':6, 'Cl':7,'Br':8, 'I':9,
               'CHI_TETRAHEDRAL_CW':11, 'CHI_TETRAHEDRAL_CCW':12,
               'SINGLE':14, 'DOUBLE':15, 'TRIPLE':16, 'AROMATIC':17,
               'STEREOE':18, 'STEREOZ':19}

    for atom in mol.GetAtoms():
        i = atom.GetIdx()

        symbol = atom.GetSymbol()
        if symbol in nf_dict:
            nfm[i, nf_dict[symbol]] = 1
        else:
            nfm[i, 10] = 1

        chiral = str(atom.GetChiralTag())
        if chiral in nf_dict:
            nfm[i, nf_dict[chiral]] = 1

        nfm[i, 13] = atom.GetFormalCharge()

        for bond in atom.GetBonds():

            stereo = str(bond.GetStereo())
            if stereo in nf_dict:
                nfm[i, nf_dict[stereo]] = 1

            bond = str(bond.GetBondType())
            if bond in nf_dict:
                nfm[i, nf_dict[bond]] += 1
                
    # Edge feature matrix
    efm = Chem.GetAdjacencyMatrix(mol).astype('float64')
    np.fill_diagonal(efm, 1)

    degmat = np.diag(np.sum(efm, axis=1))
    degmat_m0p5 = fractional_matrix_power(degmat, -0.5)
    efm = np.matmul(np.matmul(degmat_m0p5, efm), degmat_m0p5)

    n_heavyatom = len(efm)
    padding = maxatom - n_heavyatom
    efm = np.pad(efm, ((0, padding), (0, padding)), 'constant', constant_values=0.0)

    return nfm, efm

def gcgcn_input(smiles, maxatom=25):
    mol = Chem.MolFromSmiles(smiles)
    
    list_smarts = ['[CX4H3]', '[CX3H2v4]', '[CX2H1v4]', '[!R;CX4H2]',
                   '[!R;CX4H]', '[!R;CX4H0]', '[!R;CX3H1v4]', '[!R;CX3H0v4]',
                   '[!R;CX2H0;$([CX2H0](=*)(=*))]', '[!R;CX2H0;$([CX2H0](#*))]',
                   '[R;CX4H2]', '[R;CX4H]', '[R;CX4H0]', '[R;CX3H1v4]',
                   '[R;CX3H0v4]', '[R;CX2H0;$([CX2H0](=*)(=*))]',
                   '[R;CX2H0;$([CX2H0](#*))]', '[R;cX4h2]', '[R;cX4h]', '[R;cX3h1v4]',
                   '[R;cX3h0v4]', '[FX1H0]', '[ClX1H0]', '[BrX1H0]',
                   '[IX1H0]', '[OX2H1]', '[!R;OX2H0]', '[R;OX2H0]',
                   '[R;oX2h0]', '[OX1H0v2]', '[NX3H2v3]', '[NX3H1v3;!R]',
                   '[NX3H1v3;R]', '[nX3h1v3;R]', '[NX3H0v3;!R]', '[NX3H0v3;R]',
                   '[nX3h0v3;R]', '[NX2H0v3;!R]', '[NX2H0v3;R]', '[nX2h0v3;R]',
                   '[NX1H0v3]', '[SX2H1v2]', '[SX2H0v2;!R]', '[SX2H0v2;R]',
                   '[sX2h0v2;R]', '[SX1H0v2]', '[OX1H0v1]', '[NX3H0v4]']

    # Node Feature Matrix
    nfm = np.zeros((maxatom, len(list_smarts)))
    for smarts_index, smarts in enumerate(list_smarts):
        frag = Chem.MolFromSmarts(smarts)
        matched_atoms_indexes = mol.GetSubstructMatches(frag)
        nfm[matched_atoms_indexes, smarts_index] = 1

    # Edge Feature Matrix
    efm = Chem.GetAdjacencyMatrix(mol).astype('float64')
    np.fill_diagonal(efm, 1)

    degmat = np.diag(np.sum(efm, axis=1))
    degmat_m0p5 = fractional_matrix_power(degmat, -0.5)
    efm = np.matmul(np.matmul(degmat_m0p5, efm), degmat_m0p5)

    n_heavyatom = len(efm)
    padding = maxatom - n_heavyatom
    efm = np.pad(efm, ((0, padding), (0, padding)), 'constant', constant_values=0.0)

    return nfm, efm

def mol_frag(smiles, method, discriptor=0, duplicate=False):
    # Procedures are cited from https://doi.org/10.1186/s13321-019-0382-3

    # Data input
    df_smarts = pd.read_excel(here + "\\model\\fragmentation.xlsx",
                              sheet_name=method, index_col=0)
    
    # Matching all of fragments regardless of atom duplication
    mol = Chem.MolFromSmiles(smiles)
    
    matching = []
    for smarts in df_smarts['SMARTS']:
        patt = Chem.MolFromSmarts(smarts)
        matching.append( mol.GetSubstructMatches(patt) )
    
    matching = pd.Series(data=matching, index=df_smarts.index, name=smiles)
    matching = matching[matching != ()]
    
    # Sorting fragments by discriptors
    if discriptor != 0:
        discript = df_smarts.loc[matching.index].iloc[:,-discriptor:]
        discript = discript.sort_values(by=discript.columns.tolist(),
                                        axis=0, ascending=False).index.tolist()
        
        matching = matching.loc[discript]
    
    # If zero fragment found
    if len(matching) == 0:
        matching_result = pd.Series(0, index=df_smarts.index, name=smiles)
        return  matching_result
    
    # Finding fragments allowing atom duplication but not subsets
    if duplicate == True:
        flatten_matching = list(chain(*matching))
        flatten_matching = sorted(flatten_matching, key=len)
        
        atom_trashcan = []
        
        for i, candidate_atoms in enumerate(flatten_matching):
            removable_atoms = candidate_atoms
            candidate_atoms = set(candidate_atoms)
            
            for other_atoms in flatten_matching[i+1:]:
                other_atoms = set(other_atoms)
                
                if candidate_atoms != other_atoms and \
                candidate_atoms.intersection(other_atoms) == candidate_atoms:
                    atom_trashcan.append(removable_atoms)
                    break
        
        for trash_atoms in atom_trashcan:
            flatten_matching.remove(trash_atoms)
        
        matching_result = pd.Series(0, index=df_smarts.index, name=smiles)
        for smarts_n in matching.index:
            for candidate_atoms in matching[smarts_n]:
                if candidate_atoms in flatten_matching:
                    matching_result[smarts_n] += 1
        
        return matching_result
    
    # Finding fragments without atom duplication
    n_atom = mol.GetNumAtoms()
    n_frag = len( list(chain(*matching.values)) ) - 1
    
    for jump in range(n_frag):
        atom_box = ()
        matching_result = pd.Series(0, index=df_smarts.index, name=smiles)
        
        for smarts_n in matching.index:
            for candidate_atoms in matching[smarts_n]:
                if jump > 0:
                    jump -= 1
                    continue
                
                if set(atom_box).intersection(candidate_atoms) == set():
                    atom_box += candidate_atoms
                    matching_result[smarts_n] += 1
        
        atom_box = list(atom_box)
        atom_box.sort()
        
        # If fragmentation successed
        if atom_box == list(range(n_atom)):
            return matching_result
    
    # If no solution found
    matching_result = pd.Series(0, index=df_smarts.index, name=smiles)
    return matching_result
