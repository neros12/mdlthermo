import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdchem, rdMolDescriptors
from rdkit.Chem.rdchem import Mol


def Init(mol_input):
    num_heavy_atoms = Mol.GetNumHeavyAtoms(mol_input)
    dic_params = {
        "strClassL1": str(),
        "strClassL2": str(),
        "strClassL3": str(),
        "iClassL1": int(),
        "iClassL2": int(),
        "iClassL3": int(),
        "strError": str(),
        "strAtomclass": str(),
        "strMainclass": str(),
        "strSubclass ": str(),
        "iAtomclass": int(),
        "iMainclass": int(),
        "iSubclass": int(),
        # Variables
        "strTInChI": str(),  # InChI string
        "strTFormula": str(),  # Formula part in InChI
        "strTAtom": str(),  # Atom connection part in InChI
        "strTHydrogen": str(),  # Hydrogen atom part in InChI
        "strProton": str(),  # Proton part in InChI
        "strTCharge": str(),  # Chaged part in InChI
        "strInChI": str(),  # InChI string
        "strFormula": str(),  # Formula part in InChI
        "strAtom": str(),  # Atom connection part in InChI
        "strHydrogen": str(),  # Hydrogen atom part in InChI
        "strCharge": str(),  # Chaged part in InChI
        "strBond": str(),  # Bonding part In InChI
        "strTetra": str(),  # Tetra stereo part in InChI
        "strStereo": str(),  # Stereochemistry part in InChI
        "strIsotrope": str(),  # Isotrope pate in in InChI
        "strfeature": str(),  # Feature of chemical compound.
        "iNumSubStruc": int(),
        "ithSub": int(),
        "strSubInChI": [""] * 20,
        "strSubFormula": [""] * 20,
        "strSubConnect": [""] * 20,
        "strSubHydrogen": [""] * 20,
        "strSubCharge": [""] * 20,
        # Numbers of important atoms
        "iSubMetal": np.zeros(20, dtype=int),
        "iNumSubCarbon": np.zeros(20, dtype=int),
        "iNumSubNitrogen": np.zeros(20, dtype=int),
        "iNumSubOxygen": np.zeros(20, dtype=int),
        "iNumSubSilicon": np.zeros(20, dtype=int),
        "iNumSubPhosphine": np.zeros(20, dtype=int),
        "iNumSubSulfur": np.zeros(20, dtype=int),
        "iNumSubHydrogen": np.zeros(20, dtype=int),
        "iNumSubHalogen": np.zeros(20, dtype=int),  # Number of Halogen
        "iNumSubFluorine": np.zeros(20, dtype=int),  # Number of Fluorine
        "iNumSubChlorine": np.zeros(20, dtype=int),  # Number of Chlorine
        "iNumSubBromine": np.zeros(20, dtype=int),  # Number of Bromine
        "iNumSubIodine": np.zeros(20, dtype=int),  # Number of Iodine
        "iSubCountAtom": np.zeros(20, dtype=int),
        "iSubAtom": np.zeros((20, 255), dtype=int),
        "iOrgSub": np.zeros(20),
        # variables for substructure
        "iKindKGuest": int(),
        "iAtoms": np.zeros(255, dtype=int),
        "iNumRing": int(),  # Number of rings
        "iNumUnsatbond": int(),  # Number of unsaturated bonds
        "iNumTriple": int(),
        "iNumDouble": int(),
        "iNumBranch": int(),  # Number of branches
        "iNumAromatic": int(),  # Number of aromatics
        "ring": np.zeros(255, dtype=int),
        "ringatom": np.zeros((255, 255), dtype=int),
        "iringNUM": np.zeros(255, dtype=int),
        "iloc": np.zeros(20, dtype=int),
        "iBond": np.zeros(20, dtype=int),
        "iUnsatbond": np.zeros(255, dtype=int),
        "iNumbond": np.zeros(20, dtype=int),
        "iOrCheck": np.zeros(10, dtype=int),
        "FGs": np.zeros(255, dtype=int),
        "iMulSub": np.zeros(20, dtype=int),
        "nstdBond": int(),
        "iSubAnion": int(),
        "iGuestAtom": int(),
        "iCountAtom": int(),  # Total atom number in molecule
        "iKindAtom": int(),  # Total kind of atom in molecule
        "iPartNum": int(),
        "iNumCarbon": int(),  # Number of Carbon
        "iNumOxygen": int(),  # Number of Oxygen
        "iNumNitrogen": int(),  # Number of Nitrogen
        "iNumSulfur": int(),  # Number of Sulfur
        "iNumSilicon": int(),  # Number of Silicon
        "iNumPhosphine": int(),  # Number of Phosphorus
        "iNumHalogen": int(),  # Number of Halogen
        "iNumFluorine": int(),  # Number of Fluorine
        "iNumChlorine": int(),  # Number of Chlorine
        "iNumBromine": int(),  # Number of Bromine
        "iNumIodine": int(),  # Number of Iodine
        "iNumHydrogen": int(),  # Number of Hydrogen
    }

    return dic_params


if __name__ == "__main__":
    mol = Chem.MolFromSmiles("CCCC")
    mol = Chem.inchi.MolFromInchi(
        "InChI=1S/H4N2.2NO3.Ni/c1-2;2*2-1(3)4;/h1-2H2;;;/q;2*-1;+2"
    )
    dic_params = Init(mol)

    dic_params["strTFormula"] = rdMolDescriptors.CalcMolFormula(mol)


def Cal_Parts(mol, dic_params):
    return


def Do_Classification(mol):
    parameters = Init(mol)

    return
