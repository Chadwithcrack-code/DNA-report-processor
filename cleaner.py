import re

def detect_type(sequence_series):
    # Sample the first few sequences and combine to check characters safely
    sample = sequence_series.dropna().iloc[:100].str.upper().str.cat()
    if not sample: return "Unknown"
    
    # If it contains standard amino acids not found in DNA/RNA
    if bool(re.search(r"[EFILPQ]", sample)):
        return "Protein"
    # If it has Uracil, it's RNA
    elif 'U' in sample:
        return "RNA"
    # Default to DNA
    else:
        return "DNA"


def clean(df):
    # 1. Detect type once for the whole column to avoid row-by-row overhead
    mol_type = detect_type(df["sequence"])
    
    df["sequence"] = df["sequence"].astype("string")
    df["seq_length"] = df["sequence"].str.strip().str.len()
    df["molecule_type"] = mol_type
    
    # 2. Only run Nucleotide math if the sequence is DNA/RNA
    if mol_type in ["DNA", "RNA"]:
        df["atg"] = df["sequence"].str.count(r"ATG|AUG", flags=re.IGNORECASE)
        df["stop codons"] = df["sequence"].str.count(r"TAA|TAG|TGA|UAA|UAG|UGA", flags=re.IGNORECASE)
        
        # Prevent division by zero errors by replacing 0 with 1 for the denominator
        valid_len = df["seq_length"].replace(0, 1)
        valid_atg = df["atg"].replace(0, 1)
        
        df["atg_to_length(%)"] = (df["atg"] / valid_len * 100).round(3)
        df["stop_to_atg_ratio"] = (df["stop codons"] / valid_atg).round(2)
        
        gc_count = df["sequence"].str.count('G', flags=re.IGNORECASE) + df["sequence"].str.count('C', flags=re.IGNORECASE)
        df["GC (%)"] = ((gc_count / valid_len) * 100).round(3)
        df["N (%)"] = ((df["sequence"].str.count("N", flags=re.IGNORECASE) / valid_len) * 100).round(3)
    else:
        # If it's a Protein, fill with NA to keep the dataframe structure clean but save memory
        for col in ["atg", "stop codons", "atg_to_length(%)", "stop_to_atg_ratio", "GC (%)", "N (%)"]:
            df[col] = None

    return df