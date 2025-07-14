# offtarget_selectivity/01_nebula/generate_library.py

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
import os


def generate_library(input_smiles: str, num_variants: int = 20):
    print(f"Starting NEBULA generative layer for input: {input_smiles}")

    mol = Chem.MolFromSmiles(input_smiles)
    if mol is None:
        print("❌ Error: Invalid input SMILES")
        return

    mol = Chem.AddHs(mol)
    valid_mols = []

    for i in range(num_variants):
        try:
            # Generate conformer
            mol_copy = Chem.Mol(mol)
            AllChem.EmbedMolecule(mol_copy, randomSeed=42 + i)
            AllChem.UFFOptimizeMolecule(mol_copy)
            valid_mols.append(mol_copy)
        except Exception as e:
            print(f"⚠️ Skipping variant {i} due to error: {e}")

    if not valid_mols:
        print("❌ No valid molecules were generated.")
        return

    # Write valid molecules to SDF
    output_path = "nebula/generated_library.sdf"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    writer = Chem.SDWriter(output_path)

    for mol in valid_mols:
        Chem.rdDepictor.Compute2DCoords(mol)
        writer.write(mol)
    writer.close()

    print(f"✅ Generated {len(valid_mols)} 3D molecules and saved to {output_path}")
