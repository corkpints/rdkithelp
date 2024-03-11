import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd 

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Draw.MolToImage(mol)
    else:
        return None

def main():
      # Custom CSS
    st.markdown(
        f"""
        <style>
            .reportview-container .main .block-container{{
                max-width: 900px;
                padding-top: 2rem;
                padding-right: 2rem;
                padding-left: 2rem;
                padding-bottom: 2rem;
            }}
            .css-1aumxhk{{
                background-color: #FF3399;
            }}
            .css-9ck3ik{{
                color: #203864;
            }}
        </style>
        """,
        unsafe_allow_html=True
    )
    
    st.title('SMILES to Molecule Image Converter')
    smiles_input = st.text_input("Enter SMILES string:")
    df = pd.read_csv("for_app.csv")
    st.dataframe(df) 

if __name__ == "__main__":
    main()
