# utils.py
import pandas as pd
import requests
import time
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, QED, inchi
from urllib.parse import quote
from io import BytesIO

def validate_smiles(smiles):
    """Validate SMILES string using RDKit"""
    if not smiles:
        return False
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def fetch_smiles_pubchem(inchikey):
    """Fetch SMILES from PubChem using InChIKey with multiple attempt strategies"""
    try:
        # Strategy 1: Direct property lookup
        url1 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/IsomericSMILES/JSON"
        response = requests.get(url1)
        if response.status_code == 200:
            data = response.json()
            return data['PropertyTable']['Properties'][0]['IsomericSMILES']
        
        # Strategy 2: Get CID first, then SMILES
        url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
        response = requests.get(url2)
        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data:
                cid = data['IdentifierList']['CID'][0]
                smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON"
                smiles_response = requests.get(smiles_url)
                if smiles_response.status_code == 200:
                    smiles_data = smiles_response.json()
                    return smiles_data['PropertyTable']['Properties'][0]['IsomericSMILES']
        return None
    except Exception as e:
        print(f"Error fetching SMILES for {inchikey}: {str(e)}")
        return None

def generate_inchikey_rdkit(smiles):
    """Generate InChIKey from SMILES using RDKit"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            inchi_str = inchi.MolToInchi(mol)
            inchikey = inchi.InchiToInchiKey(inchi_str)
            return inchikey
        return None
    except Exception as e:
        print(f"Error generating InChIKey for SMILES {smiles}: {str(e)}")
        return None

def calculate_properties(smiles):
    """Calculate molecular properties using RDKit"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        mw = float(Descriptors.ExactMolWt(mol))
        tpsa = float(Descriptors.TPSA(mol))
        
        properties = {
            'MW': mw,
            'QED': float(QED.default(mol)),
            'HBD': int(rdMolDescriptors.CalcNumHBD(mol)),
            'HBA': int(rdMolDescriptors.CalcNumHBA(mol)),
            'Heavy Atoms': int(mol.GetNumHeavyAtoms()),
            'TPSA': tpsa,
            'PSA/MW': tpsa/mw if mw > 0 else None
        }
        
        properties['NPOL'] = properties['HBD'] + properties['HBA']
        
        return properties
    except Exception as e:
        print(f"Error calculating properties for SMILES {smiles}: {str(e)}")
        return None

def convert_df_to_csv(df):
    """Convert DataFrame to CSV"""
    return df.to_csv(index=False).encode('utf-8')

def convert_df_to_excel(df):
    """Convert DataFrame to Excel"""
    output = BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False)
    return output.getvalue()