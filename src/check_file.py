import os
files = os.listdir("../data/raw/DMS_ProteinGym_substitutions/")
# If your folder name is slightly different, use that path instead:
# files = os.listdir("./data/raw/DMS_ProteinGym_substitutions/")

print("Files available in your directory:")
for f in files[:10]: # Print the first 10 files
    print(f"- {f}")