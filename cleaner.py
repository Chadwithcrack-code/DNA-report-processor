
import re

def detect_type(df):
    # Convert a sample of the sequence to uppercase to be safe
    if df["sequence"][:1000].str.contains(r"P|V|L|I|M|F|Y|W|S|Q|K|R|H|D|E|O", flags=re.IGNORECASE, na=False).any():
        return "Protein"
    
    # Default to DNA/cDNA
    elif df["sequence"][:1000].str.contains(r"T", flags=re.IGNORECASE, na=False).any():
        return "DNA"
    
    # If it has Uracil and no Thymine, it's raw RNA
    elif df["sequence"][:1000].str.contains(r"U", flags=re.IGNORECASE, na=False).any():
        return "RNA"

    else: return "Error"


def clean(df):
    print("Ticker")
    #ATG count
    df["atg"] = df["sequence"].str.count("ATG", flags=re.IGNORECASE)
    print("atg count complete")

    #sequence length
    df["seq_length"] = df["sequence"].str.count(r"\w")
    print("sequence count complete")

    #atg to sequence length ratio
    df["atg_to_length(%)"] = df["atg"] / df["seq_length"] * 100
    df["atg_to_length(%)"] = df["atg_to_length(%)"].round(3)
    print("atg to length complete")

    #stop codon count
    df["stop codons"] = df["sequence"].str.count(r"TAA|TAG|TGA",flags=re.IGNORECASE)
    print("stop codons count complete")

    #start to stop codon ratio
    df["atg_to_stop codon"] = "1: " + ((df['stop codons']/ df["atg"]).round(2)).astype(str)
    print("atg to stop codon count complete")

    #sequence count
    df["GC (%)"] = ((df["sequence"].str.count('G', flags=re.IGNORECASE) + df["sequence"].str.count('C', flags=re.IGNORECASE)) / df["seq_length"]) * 100
    df["GC (%)"] = df["GC (%)"].round(3)
    print("GC count complete")

    #N to sequence ratio
    df["N (%)"] = (df["sequence"].str.count("N", flags=re.IGNORECASE) / df["seq_length"]) * 100
    df["N (%)"] = df["N (%)"].round(3)
    print("N count complete")

    return df

