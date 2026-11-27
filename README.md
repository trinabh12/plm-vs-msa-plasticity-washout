# PLM vs MSA: Plasticity and Washout in Protein Variant Analysis

A comparative research repository for studying protein variant effects with **protein language models (PLMs)** and **multiple-sequence-alignment (MSA)-based representations**. The project combines deep mutational scanning (DMS) data, locally extracted ESM-2 embeddings, protein structure prediction, structural confidence metrics, and visualization workflows.

The supplied example focuses on **BRCA1 (UniProt P38398)** and the BRCA1 assay from `ProteinGym_Findlay_2018`, but the data layout is designed to support many ProteinGym targets.

## Research Questions

- How do PLM representations encode residue-level and variant-level protein information?
- How well do PLM and MSA-derived signals capture experimentally measured variant effects?
- Do protein variants alter predicted structure, confidence, or local geometry?
- Can structural and embedding-based analyses help explain stability or "plasticity/washout" effects?

## Analysis Pipeline

```text
ProteinGym DMS CSV + MSA weights
              |
              v
Mutant sequence preparation
              |
      +-------+-------+
      |               |
      v               v
ESM-2 embeddings   Structure prediction
(650M, local)      (ESMFold / AlphaFold)
      |               |
      v               v
PCA / t-SNE / UMAP  pLDDT, RMSD, Ramachandran
      |               |
      +-------+-------+
              |
              v
Stability, annotation, and explainability analysis
```

## Repository Layout

```text
.
├── data/
│   ├── raw/
│   │   ├── DMS_ProteinGym_substitutions/  # ProteinGym DMS assay CSVs
│   │   ├── DMS_msa_weights/               # Precomputed MSA-derived NPY weights
│   │   └── ProteinGym_AF2_structures/     # AlphaFold structures bundled for targets
│   └── processed/
│       ├── embeddings/                    # ESM-2 embedding matrices (.npy)
│       └── predicted_structures/
│           ├── esmfold/                   # ESMFold predictions and alignments
│           └── alphafold/                 # Downloaded AlphaFold entries
├── notebooks/
│   ├── structure_prediction.ipynb
│   ├── structural_annotation.ipynb
│   ├── stability_assessment.ipynb
│   ├── benchmarking_and_embeddings.ipynb
│   └── explainability_and_stability.ipynb
├── src/
│   ├── data_loader.py                     # DMS assay loading and variant selection
│   ├── extract_embeddings.py              # ESM-2 local embedding extraction
│   ├── structure_engine.py                # AlphaFold/ESMFold fetching and structure metrics
│   └── check_file.py                      # Utility for inspecting available assay files
├── visualisations/
│   ├── final_superposition.png
│   └── render_figure.pml                  # PyMOL rendering script
├── environment.yml                        # Reproducible Conda environment
└── README.md
```

## Data Overview

### Deep mutational scanning data

ProteinGym substitution assay CSVs are stored in `data/raw/DMS_ProteinGym_substitutions/`. A typical file contains:

| Column | Description |
|---|---|
| `mutant` | Variant identifier, for example `M1I` |
| `mutated_sequence` | Full amino-acid sequence containing the variant |
| `DMS_score` | Experimental assay or fitness score |
| `DMS_score_bin` | Optional binary or thresholded assay label |

The repository includes 217 such assay files. The loader expects `mutant` and `mutated_sequence` columns.

### MSA weights

The `data/raw/DMS_msa_weights/` directory contains precomputed `.npy` files used as MSA/indirect evolutionary baselines or features. These arrays should be loaded and interpreted consistently with the target assay and ProteinGym conventions.

### Structures

- `data/raw/ProteinGym_AF2_structures/` contains reference AlphaFold PDB files.
- `data/processed/predicted_structures/alphafold/` contains dynamically downloaded AlphaFold database entries.
- `data/processed/predicted_structures/esmfold/` contains ESMFold predictions and superposition-aligned PDB files.

## Environment Setup

### Requirements

- Windows, Linux, or macOS
- [Conda](https://docs.conda.org/) or [Miniconda](https://docs.anaconda.com/miniconda/)
- Internet access for Hugging Face model downloads and structure API requests
- An NVIDIA GPU is recommended for ESM-2 embedding extraction; CPU execution is supported but slower

### Create the environment

From the repository root:

```powershell
conda env create -f environment.yml
conda activate protein_ai
python -m ipykernel install --user --name protein_ai --display-name "Python (protein_ai)"
```

If the environment already exists, update it with:

```powershell
conda env update -f environment.yml --prune
```

Key dependencies include PyTorch, Transformers, Hugging Face Accelerate, pandas, NumPy, scikit-learn, Matplotlib, Seaborn, UMAP, Biopython, and PyMOL Open Source.

## Running the Workflow

### 1. Use the Jupyter notebooks

The easiest entry point is the notebook sequence:

```powershell
jupyter lab
```

Open the notebooks in this order:

1. `structure_prediction.ipynb` — load assay variants, fetch/predict structures, calculate RMSD, produce geometry plots, and align structures.
2. `structural_annotation.ipynb` — retrieve UniProt annotations and map domains or functional regions onto a protein sequence.
3. `benchmarking_and_embeddings.ipynb` — pool ESM-2 embeddings and compare variants through PCA, t-SNE, and UMAP.
4. `explainability_and_stability.ipynb` — inspect ESM attention and mutation-effect proxies.
5. `stability_assessment.ipynb` — compile stability-change results into a standardized report.

Launch Jupyter from the **repository root**. The notebook paths assume the root-relative `../src` and `../data` layout.

### 2. Extract ESM-2 embeddings locally

Change into `src/` because the script uses paths relative to that directory:

```powershell
cd src
python extract_embeddings.py
```

The script loads `facebook/esm2_t33_650M_UR50D`, truncates sequences longer than 1022 residues, removes CLS/EOS tokens, and saves a residue-by-1280 hidden-state matrix for each variant:

```text
data/processed/embeddings/<variant_id>_embeddings.npy
```

By default, the script uses `TARGET_ASSAY = "BRCA1_HUMAN"` and processes ten variants. Update that constant in `src/extract_embeddings.py` to process another assay prefix.

### 3. Fetch or predict structures

Example use from the `src/` directory:

```python
from data_loader import load_dms_assay
from structure_engine import predict_esmfold, fetch_alphafold_db

variants = load_dms_assay("BRCA1_HUMAN", max_mutants=5)
fetch_alphafold_db("P38398")

for variant_id, sequence in variants.items():
    predict_esmfold(variant_id, sequence)
```

The ESMFold API is limited to the first 150 residues in `predict_esmfold()` to reduce server overload and timeout risk. Existing PDB outputs are reused as a local cache.

## Core Python Modules

### `src/data_loader.py`

- Finds ProteinGym CSV files by assay name.
- Requires `mutant` and `mutated_sequence` columns.
- Returns a dictionary of variant IDs to mutated sequences.
- Limits variant count through `max_mutants` for exploratory runs.

### `src/structure_engine.py`

- Fetches AlphaFold PDB structures through the AlphaFold API.
- Queries the ESM Atlas ESMFold API with retries for HTTP 504 responses.
- Extracts C-alpha pLDDT values from PDB B-factors.
- Aligns structures and calculates C-alpha RMSD with Biopython.

### `src/extract_embeddings.py`

- Loads ESM-2 650M and uses GPU automatically when available.
- Produces a residue-by-hidden-state NumPy matrix for each variant.
- Saves processed embeddings under `data/processed/embeddings/`.

## Output and Visualization

- **pLDDT** summarizes local predicted structure confidence.
- **C-alpha RMSD** measures structural deviation after superposition.
- **Ramachandran plots** visualize backbone geometry.
- **PCA, t-SNE, and UMAP** visualize relationships among pooled embedding vectors.
- **PyMOL** and `py3Dmol` support static or interactive structural comparison.
- `visualisations/final_superposition.png` provides an example final structure superposition.
- `visualisations/render_figure.pml` contains a PyMOL figure-rendering recipe.

## Reproducibility Notes

- Pin analysis inputs by target assay, variant subset, and random seed.
- Keep raw DMS CSVs and downloaded structures unchanged; place derived outputs in `data/processed/`.
- Record the Python/PyTorch/CUDA versions used for GPU runs.
- Transformer weights and API-predicted structures can change over time; cache and version outputs for repeatable comparisons.
- For the sample BRCA1 workflow, verify that assay filenames and UniProt IDs match the target being analyzed.
- Relative paths in the Python modules assume the `src/` directory is the current working directory for scripts and the repository root is the working directory for notebooks.

## Troubleshooting

### CUDA is not detected

```powershell
python -c "import torch; print(torch.cuda.is_available()); print(torch.__version__)"
```

Install/reinstall a PyTorch build matching the available NVIDIA driver and CUDA toolkit. CPU execution remains available.

### Hugging Face model download fails

- Verify internet access and Hugging Face connectivity.
- Authenticate if the model or network requires credentials.
- Ensure enough disk space is available for the ESM-2 650M checkpoint.

### ESMFold or AlphaFold API requests fail

- Re-run structure fetching; ESMFold already retries HTTP 504 responses.
- Check whether a PDB cache file already exists.
- The external APIs may be unavailable or rate-limited.

### A DMS assay is not found

List the available files or adjust the assay prefix:

```powershell
cd src
python check_file.py
```

Use the matching part of a filename before the `.csv` extension, for example `BRCA1_HUMAN` for `BRCA1_HUMAN_Findlay_2018.csv`.

### A notebook cannot import `src` modules

Start Jupyter from the repository root and run the notebook setup cell that adds `../src` to `sys.path`. If using VS Code, select the `protein_ai` kernel.

## Citation and Data Sources

Please cite the original ProteinGym, DMS, UniProt, AlphaFold, ESMFold, ESM-2, and assay publications that correspond to the data or method used in a downstream analysis.

This repository does not include an explicit license file. Check the terms of the bundled datasets and original source services before redistribution or publication.

## Contributing

Contributions should preserve the raw/processed data split, keep large derived artifacts out of source-control patches where possible, and document any new assay naming conventions or notebook outputs.
