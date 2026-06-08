
from Bio.Seq import Seq
from cleaner import detect_type
import pandas as pd

def analyze_amino(df):
    #percentage of unusual atg to sequence length (this can be changed with your discretion)
    percent = 3
    concentrated_atg = []

    #pcr = Protein coding regions: this is a df that only includes the protein coding regions
    pcr = df[((df["stop codons"] / df["atg"]) < 2) & (df["atg_to_length(%)"] > percent)]

    #Raw DNA strand for processing
    raw_seq = pcr.loc[:, "sequence"]

    #Converting the raw DNA strandinto a Biopython Seq object so that we use .translate() method
    #The translate method adds a * by the stop codon which we can use split to get the protein
    #trimmed the sequence to be a multiple of 3 (codon length)

    type = detect_type(df)
    print(type)
    try:
        if type == "DNA":
            for x in raw_seq:
                concentrated_atg.extend(str(Seq(x).translate()).split("*"))
        elif type == "RNA":
            for x in raw_seq:
                concentrated_atg.extend(str(Seq(x).back_transcribe().translate()).split("*"))
        elif type == "Protein":
            for x in raw_seq:
                concentrated_atg.extend(str(x).split("*"))
        else: raise ValueError
    except ValueError:
        return "problem with format type:"
    




    #Now we filter out for actual protein, which must start with Methionine(M)
    #count to label the proteins
    count = 0
    c_genes = pd.Series([])
    for x in concentrated_atg:
        if x.startswith("M"):
            c_genes[f"Proteins_{count}"] = x
            count += 1
            print(x)
    #Next we take all the proteins and join them into a coherent string 
    #We break up element into a python list
    #finally we use turn it into a pd.Series so that we can calculate the percetage using value_counts
    aa_percentages = pd.Series(list(c_genes.str.cat())).value_counts(normalize=True) * 100

    #Cleaning the DataFrame
    final_composition_df = aa_percentages.reset_index()
    final_composition_df.index += 1
    final_composition_df.columns = ["Amino Acid", "Concentration (%)"]

    # Round the numbers so it looks clean for the README screenshot
    final_composition_df["Concentration (%)"] = final_composition_df["Concentration (%)"].round(2)
    return final_composition_df

    print("AMINO ACID COMPOSITION:")
    final_composition_df.head(15)
