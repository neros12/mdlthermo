import os
from os.path import join as opj

import numpy as np
import pandas as pd
from rdkit import Chem
from pathlib import Path


df = pd.read_csv(os.path.join(Path(__file__).parents[0], "components_list.csv"))
InChIKey_to_COSMO = dict(zip(df["InChIKey"], df["Index"]))


def get_distance(coordinate1, cooredinate2):

    return np.linalg.norm(np.array(coordinate1) - np.array(cooredinate2))


def get_COSMO_file_name(SMILES: str) -> str:
    mol = Chem.MolFromSmiles(SMILES)
    InChIKey = Chem.inchi.MolToInchiKey(mol)
    COSMO_file_name = InChIKey_to_COSMO[InChIKey]

    return COSMO_file_name


def read_COSMO(SMILES: str):
    A = 0.0  # Total surface area of cavity (A**2)
    V = 0.0  # Total volume of cavity (A**3)
    atoms = []
    segments = []

    COSMO_file_name = get_COSMO_file_name(SMILES)
    COSMO_file_path = opj(
        Path(__file__).parent,
        "cosmo_files",
        f"{COSMO_file_name}.cosmo",
    )

    with open(COSMO_file_path, "r") as cosmo_file:
        isCoordinate = False
        isSegment = False
        for line in cosmo_file.readlines():
            if "Total surface area of cavity (A**2)" in line:
                A = float(line.split("=")[1].strip())

            elif "Total volume of cavity (A**3)" in line:
                V = float(line.split("=")[1].strip())

            elif "!DATE" in line:
                isCoordinate = True
                continue

            elif isCoordinate:
                if "end" in line:
                    isCoordinate = False
                    continue

                (
                    _,
                    x_coordinate,
                    y_coordinate,
                    z_coordinate,
                    _,
                    _,
                    _,
                    atom_type,
                    _,
                ) = line.split()
                atoms.append(
                    {
                        "atom_index": len(atoms),
                        "atom_type": atom_type,
                        "coordinate": {
                            "x": float(x_coordinate),
                            "y": float(y_coordinate),
                            "z": float(z_coordinate),
                        },
                    }
                )

            elif "position (X, Y, Z) [au]" in line:
                isSegment = True
                continue

            elif isSegment:
                if "\n" == line:
                    continue

                (
                    segment_index,
                    atom_index,
                    x_position,
                    y_position,
                    z_position,
                    charge,
                    area,
                    _,
                    potential,
                ) = line.split()
                segments.append(
                    {
                        "segment_index": int(segment_index) - 1,
                        "atom_index": int(atom_index) - 1,
                        "coordinate": {
                            "x": float(x_position),
                            "y": float(y_position),
                            "z": float(z_position),
                        },
                        "charge": float(charge),
                        "area": float(area),
                        "potential": float(potential),
                    }
                )

    return A, V, atoms, segments
