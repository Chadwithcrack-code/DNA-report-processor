
from Bio.Seq import Seq
from cleaner import detect_type
import pandas as pd

def analyze_amino(df, perc, min_amino_length=50, tb=1):
    #percentage of unusual atg to sequence length (this can be changed with your discretion)
    concentrated_atg = []

    #pcr = Protein coding regions: this is a df that only includes the protein coding regions
    pcr = df[((df["stop codons"] / df["atg"]) < 2) & (df["atg_to_length(%)"] > perc)]

    #Raw DNA strand for processing
    raw_seq = pcr.loc[:, "sequence"]

    #Converting the raw DNA strandinto a Biopython Seq object so that we use .translate() method
    #The translate method adds a * by the stop codon which we can use split to get the protein
    #trimmed the sequence to be a multiple of 3 (codon length)

    #type = detect_type(df)
    #print(type)
    # Initialize a list to hold the mapped data instead of a flat list
    extracted_proteins = []
    count = 0

    molecule_type = detect_type(df)

    # .items() yields BOTH the original pandas index and the sequence string
    for original_index, seq_value in raw_seq.items():
        
        # 1. Translate and split based on molecule type
        if molecule_type == "DNA":
            chunks = str(Seq(seq_value).translate()).split("*")
        elif molecule_type == "RNA":
            chunks = str(Seq(seq_value).back_transcribe().translate(table=tb)).split("*")
        elif molecule_type == "Protein":
            chunks = str(seq_value).split("*")
        else: 
            raise ValueError("problem with format type:")

        # 2. Immediately filter for Methionine while we still know the original_index
        for chunk in chunks:
            if chunk.startswith("M"):
                amino_length = len(chunk)

                if chunk.startswith("M") and amino_length >= min_amino_length:
                    extracted_proteins.append({
                        "original_index": original_index,
                        "protein_label": f"Proteins_{count}",
                        "sequence": chunk,
                        "length": amino_length
                    })
                    count += 1
                    

    # 3. Convert the structured list back into a Pandas DataFrame
    c_genes_df = pd.DataFrame(extracted_proteins)
    if c_genes_df.empty:
        return pd.DataFrame(columns=["Amino Acid", "Concentration (%)"]), c_genes_df

    # 4. Set the index of this new table to perfectly match your main Python series
    if not c_genes_df.empty:
        c_genes_df.set_index("original_index", inplace=True)
    
    print(c_genes_df)
    #Next we take all the proteins and join them into a coherent string 
    #We break up element into a python list
    #finally we use turn it into a pd.Series so that we can calculate the percetage using value_counts
    aa_percentages = pd.Series(list(c_genes_df["sequence"].str.cat())).value_counts(normalize=True) * 100

    #Cleaning the DataFrame
    final_composition_df = aa_percentages.reset_index()
    final_composition_df.index += 1
    final_composition_df.columns = ["Amino Acid", "Concentration (%)"]

    # Round the numbers so it looks clean for the README screenshot
    final_composition_df["Concentration (%)"] = final_composition_df["Concentration (%)"].round(2)
    return final_composition_df, c_genes_df
