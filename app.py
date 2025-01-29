# app.py
import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import inchi
from utils import (fetch_smiles_pubchem, generate_inchikey_rdkit, 
                  calculate_properties, convert_df_to_csv, convert_df_to_excel,
                  validate_smiles)
import time
from io import BytesIO

st.set_page_config(page_title="Chemical Structure Identifier", page_icon="🧪", layout="wide")

def initialize_session_state():
    if 'processed_df' not in st.session_state:
        st.session_state.processed_df = None
    if 'manual_inputs' not in st.session_state:
        st.session_state.manual_inputs = {}
    if 'current_process_state' not in st.session_state:
        st.session_state.current_process_state = 'initial'
    if 'missing_smiles_queue' not in st.session_state:
        st.session_state.missing_smiles_queue = []

def process_chemical_identifiers(df):
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Copy the dataframe
    if st.session_state.processed_df is None:
        processed_df = df.copy()
    else:
        processed_df = st.session_state.processed_df.copy()
    
    total_rows = len(df)
    current_row = 0
    missing_smiles = []
    
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
                missing_smiles.append({
                    'idx': idx,
                    'name': df.loc[idx, 'Name'],
                    'inchikey': inchikey
                })
        
        # Case 2: Has SMILES but no InChIKey
        elif pd.notna(smiles) and (pd.isna(inchikey) or inchikey == ''):
            status_text.text(f"Processing row {current_row}/{total_rows}: Generating InChIKey using RDKit")
            new_inchikey = generate_inchikey_rdkit(smiles)
            if new_inchikey:
                processed_df.loc[idx, 'InChIKey'] = new_inchikey
        
        time.sleep(0.1)  # Prevent rate limiting
    
    progress_bar.progress(1.0)
    status_text.text("Initial processing complete!")
    
    st.session_state.processed_df = processed_df
    st.session_state.missing_smiles_queue = missing_smiles
    st.session_state.current_process_state = 'manual_input'
    
    return processed_df, missing_smiles

def handle_manual_input(missing_smiles):
    if not missing_smiles:
        st.session_state.current_process_state = 'completed'
        return st.session_state.processed_df
    
    st.write("### Manual SMILES Input Required")
    st.write("Some compounds need manual SMILES input. Please enter the SMILES for each compound:")
    
    for item in missing_smiles:
        col1, col2, col3 = st.columns([3, 2, 1])
        
        with col1:
            st.write(f"Compound: {item['name']}")
            st.write(f"InChIKey: {item['inchikey']}")
        
        with col2:
            key = f"manual_smiles_{item['idx']}"
            manual_smiles = st.text_input("Enter SMILES:", key=key)
            if manual_smiles:
                st.session_state.manual_inputs[item['idx']] = manual_smiles
        
        with col3:
            if st.button("Confirm", key=f"confirm_{item['idx']}"):
                if validate_smiles(manual_smiles):
                    st.session_state.processed_df.loc[item['idx'], 'SMILES'] = manual_smiles
                    st.success("SMILES added successfully!")
                else:
                    st.error("Invalid SMILES string. Please check and try again.")
    
    if st.button("Complete Manual Input"):
        st.session_state.current_process_state = 'completed'
        return st.session_state.processed_df
    
    return st.session_state.processed_df

# app.py (Add this function before the main() function)

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
#Main App
def main():
    st.title("Chemical Structure Identifier and Property Calculator")
    initialize_session_state()
    
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
        if st.session_state.current_process_state == 'initial':
            if st.button("Process SMILES and InChIKeys"):
                processed_df, missing_smiles = process_chemical_identifiers(df)
                if not missing_smiles:
                    st.session_state.current_process_state = 'completed'
        
        # Handle manual input if needed
        if st.session_state.current_process_state == 'manual_input':
            processed_df = handle_manual_input(st.session_state.missing_smiles_queue)
        
        # Show download buttons and property calculation option
        if st.session_state.current_process_state == 'completed':
            st.write("### Processed Data")
            st.dataframe(st.session_state.processed_df.head())
            
            col1, col2 = st.columns(2)
            with col1:
                csv = convert_df_to_csv(st.session_state.processed_df)
                st.download_button(
                    label="Download Processed Data (CSV)",
                    data=csv,
                    file_name="processed_chemicals.csv",
                    mime="text/csv"
                )
            
            with col2:
                excel = convert_df_to_excel(st.session_state.processed_df)
                st.download_button(
                    label="Download Processed Data (Excel)",
                    data=excel,
                    file_name="processed_chemicals.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
            
            if st.button("Calculate Properties"):
                df_with_properties = calculate_all_properties(st.session_state.processed_df)
                st.write("### Data with Properties")
                st.dataframe(df_with_properties.head())
                
                col3, col4 = st.columns(2)
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

if __name__ == "__main__":
    main()