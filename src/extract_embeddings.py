import os
import torch
import numpy as np
from transformers import AutoTokenizer, EsmModel
from data_loader import load_dms_assay # Import our new loader!

OUTPUT_EMBED_PATH = "../data/processed/embeddings/"
os.makedirs(OUTPUT_EMBED_PATH, exist_ok=True)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Executing pipeline locally. Hardware targeted: {device}")

# Load ESM-2 
model_name = "facebook/esm2_t33_650M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = EsmModel.from_pretrained(model_name).to(device)
model.eval()

def generate_local_features(protein_id, sequence):
    # ESM-2 max input length is 1022 amino acids
    # Longer sequences must be truncated
    if len(sequence) > 1022:
        sequence = sequence[:1022]
        
    inputs = tokenizer(sequence, return_tensors="pt", padding=False, truncation=True).to(device)
    
    with torch.no_grad():
        outputs = model(**inputs)
        
    hidden_states = outputs.last_hidden_state.squeeze(0).cpu().numpy()
    
    # Strip <cls> and <eos> tokens
    clean_matrix = hidden_states[1:-1, :]
    
    destination = os.path.join(OUTPUT_EMBED_PATH, f"{protein_id}_embeddings.npy")
    np.save(destination, clean_matrix)
    print(f"SUCCESS: Captured matrix for {protein_id} | Shape: {clean_matrix.shape}")

if __name__ == "__main__":
    print("\n--- Initializing Local Feature Extraction Pipeline ---")
    
    # 1. Look inside your 'data/raw/proteingym_substitutions/' folder.
    # 2. Pick a target assay file (e.g., 'BRCA1_HUMAN_Findlay_2018.csv')
    # 3. Enter part of that filename below:
    TARGET_ASSAY = "BRCA1_HUMAN" 
    
    try:
        # Load 10 variants from the CSV
        sequences_to_process = load_dms_assay(TARGET_ASSAY, max_mutants=10)
        
        # Run inference loop
        for target_id, raw_seq in sequences_to_process.items():
            generate_local_features(target_id, raw_seq)
            
    except Exception as e:
        print(f"Pipeline failed: {e}")
        print("Make sure you update TARGET_ASSAY to match a file in your raw data folder!")