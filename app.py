import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Draw.MolToImage(mol)
    else:
        return None

def main():
    st.title('SMILES to Molecule Image Converter')
    
    # Input field for SMILES string
    smiles_input = st.text_input("Enter SMILES string:")
    
    # Generate molecule image
    if st.button("Generate Molecule Image"):
        if smiles_input.strip() != "":
            molecule_image = draw_molecule(smiles_input)
            if molecule_image is not None:
                st.image(molecule_image, caption='Molecule', use_column_width=True)
            else:
                st.error("Invalid SMILES string. Please enter a valid SMILES string.")
        else:
            st.warning("Please enter a SMILES string.")

if __name__ == "__main__":
    main()
