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
    smiles_input = st.text_input("Enter SMILES string:")
    df = pd.read_csv("for_app.csv")
    st.dataframe(df) 

if __name__ == "__main__":
    main()
