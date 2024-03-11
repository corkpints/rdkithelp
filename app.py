import streamlit as st
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
import pandas as pd 

def calculate_similarity(input_smiles, smiles):
    input_mol = Chem.MolFromSmiles(input_smiles)
    mol = Chem.MolFromSmiles(smiles)
    if input_mol is None or mol is None:
        return None
    fp1 = MACCSkeys.GenMACCSKeys(input_mol)
    fp2 = MACCSkeys.GenMACCSKeys(mol)
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity
    
def colorize_similarity(similarity):
    if similarity is None:
        return ''
    elif similarity >= 0.5:
        # Green shade for high similarity
        return f'background-color: rgb({int(255 - 255 * similarity)}, 255, {int(255 - 255 * similarity)})'
    else:
        # Red shade for low similarity
        return f'background-color: rgb(255, {int(255 * similarity)}, {int(255 * similarity)})'
def main():
    st.title('SMILES to Molecule Image Converter')
    smiles_input = st.text_input("Enter SMILES string:")
    df = pd.read_csv("for_app.csv", index_col=0)
    
    if smiles_input.strip() != "":
        df['Similarity'] = df['Parent compound SMILES'].apply(lambda x: calculate_similarity(smiles_input, x))
        df = df.sort_values(by='Similarity', ascending=False)
        styled_df = df.style.applymap(colorize_similarity, subset=['Similarity'])
        # Display top ten molecules and their similarities
        st.write("Top 10 Similar Molecules:")
        st.dataframe(df.head(10))
    else:
        st.warning("Please enter a SMILES string.")

if __name__ == "__main__":
    main()
