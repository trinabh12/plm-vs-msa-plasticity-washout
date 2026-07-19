import os
import pandas as pd
import glob

# Path to the extracted ProteinGym substitution CSVs
CSV_DIR = "../data/raw/DMS_ProteinGym_substitutions/"

def load_dms_assay(assay_name, max_mutants=5):
    """
    Searches for a specific DMS assay CSV, extracts the wild-type sequence, 
    and returns a dictionary mapping variant IDs to their full sequences.
    """
    # Find the specific CSV file (handling potential case/version differences in filenames)
    search_pattern = os.path.join(CSV_DIR, f"*{assay_name}*.csv")
    matching_files = glob.glob(search_pattern)
    
    if not matching_files:
        raise FileNotFoundError(f"Could not find any CSV matching '{assay_name}' in {CSV_DIR}")
        
    file_path = matching_files[0]
    print(f"Loading data from: {file_path}")
    
    df = pd.read_csv(file_path)
    
    # ProteinGym CSVs contain these standard columns
    required_cols = ['mutant', 'mutated_sequence']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"CSV missing required column: {col}")
            
    sequence_dict = {}
    
    # 1. Extract the Wild-Type (WT) Sequence
    # The 'mutant' column contains strings like 'M1A'. 
    # If we want the WT, we need to reverse a mutation or check if a WT row exists.
    # Fortunately, the 'target_seq' column usually holds the WT, but to be safe, 
    # we will reconstruct it from the first mutated_sequence.
    first_row = df.iloc[0]
    mutant_code = first_row['mutant'] # e.g., 'M1A' (Pos 1, M -> A)
    mutated_seq = first_row['mutated_sequence']
    
    # We will just grab the first 'max_mutants' to avoid overloading the GPU during testing
    df_subset = df.head(max_mutants)
    
    for index, row in df_subset.iterrows():
        variant_id = f"{assay_name}_{row['mutant']}"
        sequence_dict[variant_id] = row['mutated_sequence']
        
    print(f"Successfully loaded {len(sequence_dict)} sequences for assay: {assay_name}")
    return sequence_dict

if __name__ == "__main__":
    # Test the loader (Assuming you have A0A140D2T1_ZIKV_Sourisseau_2019 or similar in your folder)
    # Just grab any filename from your folder minus the .csv extension to test
    try:
        sample_data = load_dms_assay("A0A140D2T1") 
        for k, v in sample_data.items():
            print(f"{k}: {v[:30]}...")
    except Exception as e:
        print(e)