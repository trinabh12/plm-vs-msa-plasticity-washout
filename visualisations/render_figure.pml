# render_final_figure.pml for Visualization
# Load the AlphaFold and ESMFold structures
load ../data/processed/predicted_structures/alphafold/AF_P38398.pdb, wt_alphafold
load ../data/processed/predicted_structures/esmfold/BRCA1_HUMAN_M1I_esmfold.pdb, mut_esmfold

# Superimpose the mutant onto the wild-type to calculate structural deviation
super mut_esmfold, wt_alphafold

# Apply professional publication colors and styling
color marine, wt_alphafold
color tv_orange, mut_esmfold
show cartoon
hide lines
bg_color white

# Center the camera strictly on the mutant fragment
zoom mut_esmfold, buffer=10

# Ray tracing for the final report figure
set ray_trace_mode, 1
set antialias, 2
png final_superposition.png, dpi=300, width=1920, height=1080