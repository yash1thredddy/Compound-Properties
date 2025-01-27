# app.py
import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import inchi
from utils import (fetch_smiles_pubchem, generate_inchikey_rdkit, 
                  calculate_properties, convert_df_to_csv, convert_df_to_excel)
import time
from io import BytesIO

st.set_page_config(page_title="Chemical Structure Identifier", page_icon="🧪", layout="wide")

def main():
    st.title("Chemical Structure Identifier and Property Calculator")
    st.write("Upload your file with chemical structures to process SMILES and InChIKeys")

    uploaded_file = st.file_uploader("Choose a file", type=['csv', 'xlsx'])
    
    if uploaded_file is not None:
        # Read the file
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        else:
            df = pd.read_excel(uploaded_file)
        
        st.write("Original Data Preview:")
        st.dataframe(df.head())
        
        # Process the data
        if st.button("Process SMILES and InChIKeys"):
            df = process_chemical_identifiers(df)
            
            # Create columns for download buttons
            col1, col2 = st.columns(2)
            
            # Download buttons for processed data
            with col1:
                csv = convert_df_to_csv(df)
                st.download_button(
                    label="Download Processed Data (CSV)",
                    data=csv,
                    file_name="processed_chemicals.csv",
                    mime="text/csv"
                )
            
            with col2:
                excel = convert_df_to_excel(df)
                st.download_button(
                    label="Download Processed Data (Excel)",
                    data=excel,
                    file_name="processed_chemicals.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
            
            # Calculate properties button
            if st.button("Calculate Properties"):
                df_with_properties = calculate_all_properties(df)
                
                st.write("Data with Properties Preview:")
                st.dataframe(df_with_properties.head())
                
                # Create columns for property data download buttons
                col3, col4 = st.columns(2)
                
                # Download buttons for data with properties
                with col3:
                    csv_props = convert_df_to_csv(df_with_properties)
                    st.download_button(
                        label="Download Data with Properties (CSV)",
                        data=csv_props,
                        file_name="chemicals_with_properties.csv",
                        mime="text/csv"
                    )
                
                with col4:
                    excel_props = convert_df_to_excel(df_with_properties)
                    st.download_button(
                        label="Download Data with Properties (Excel)",
                        data=excel_props,
                        file_name="chemicals_with_properties.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )

def process_chemical_identifiers(df):
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Copy the dataframe
    processed_df = df.copy()
    
    total_rows = len(df)
    current_row = 0
    
    for idx in df.index:
        current_row += 1
        progress = current_row / total_rows
        progress_bar.progress(progress)
        
        inchikey = df.loc[idx, 'InChIKey'] if 'InChIKey' in df.columns else None
        smiles = df.loc[idx, 'SMILES'] if 'SMILES' in df.columns else None
        
        # Case 1: Has InChIKey but no SMILES
        if pd.notna(inchikey) and (pd.isna(smiles) or smiles == ''):
            status_text.text(f"Processing row {current_row}/{total_rows}: Fetching SMILES from PubChem")
            new_smiles = fetch_smiles_pubchem(inchikey)
            if new_smiles:
                processed_df.loc[idx, 'SMILES'] = new_smiles
            else:
                # Ask user for manual input
                st.warning(f"Could not fetch SMILES for InChIKey: {inchikey}")
                manual_smiles = st.text_input(f"Enter SMILES manually for {df.loc[idx, 'Name']} (InChIKey: {inchikey})")
                if manual_smiles:
                    processed_df.loc[idx, 'SMILES'] = manual_smiles
        
        # Case 2: Has SMILES but no InChIKey
        elif pd.notna(smiles) and (pd.isna(inchikey) or inchikey == ''):
            status_text.text(f"Processing row {current_row}/{total_rows}: Generating InChIKey using RDKit")
            new_inchikey = generate_inchikey_rdkit(smiles)
            if new_inchikey:
                processed_df.loc[idx, 'InChIKey'] = new_inchikey
                
        # Case 3: If there's an InChI value instead of SMILES or InChIKey
        elif 'InChI' in df.columns and pd.notna(df.loc[idx, 'InChI']):
            status_text.text(f"Processing row {current_row}/{total_rows}: Converting InChI to SMILES and InChIKey")
            mol = Chem.MolFromInchi(df.loc[idx, 'InChI'])
            if mol:
                processed_df.loc[idx, 'SMILES'] = Chem.MolToSmiles(mol)
                inchi_str = Chem.MolToInchi(mol)
                processed_df.loc[idx, 'InChIKey'] = Chem.InchiToInchiKey(inchi_str)
        
        time.sleep(0.1)  # Prevent rate limiting
    
    progress_bar.progress(1.0)
    status_text.text("Processing complete!")
    
    return processed_df

def calculate_all_properties(df):
    processed_df = df.copy()
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    total_rows = len(df)
    current_row = 0
    
    # Initialize property columns
    property_columns = ['MW', 'QED', 'HBD', 'HBA', 'Heavy Atoms', 'NPOL', 'TPSA', 'PSA/MW']
    for col in property_columns:
        processed_df[col] = None
    
    for idx in df.index:
        current_row += 1
        progress = current_row / total_rows
        progress_bar.progress(progress)
        
        smiles = df.loc[idx, 'SMILES']
        if pd.notna(smiles):
            status_text.text(f"Calculating properties for compound {current_row}/{total_rows}")
            properties = calculate_properties(smiles)
            if properties:
                for prop, value in properties.items():
                    processed_df.loc[idx, prop] = value
    
    # Round numeric columns
    numeric_columns = ['MW', 'QED', 'TPSA', 'PSA/MW']
    for col in numeric_columns:
        if col in processed_df.columns:
            processed_df[col] = pd.to_numeric(processed_df[col], errors='coerce')
            if col == 'PSA/MW':
                processed_df[col] = processed_df[col].round(3)
            else:
                processed_df[col] = processed_df[col].round(2)
    
    progress_bar.progress(1.0)
    status_text.text("Property calculation complete!")
    
    return processed_df

if __name__ == "__main__":
    main()