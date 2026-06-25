from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from collections import Counter
import pandas as pd


def analyze_amino(df, perc, min_amino_length=50, tb=1):
    if df.empty:
        return pd.DataFrame(columns=["Amino Acid", "Concentration (%)"]), pd.DataFrame()

    molecule_type = df["molecule_type"].iloc[0] if "molecule_type" in df.columns else "DNA"
    extracted_proteins = []
    count = 0
    
    # 1. Filter Logic based on Molecule Type
    if molecule_type in ["DNA", "RNA"]:
        # This allows large genomes (which naturally have a ~3.0 ratio) to pass to the translation engine.
        pcr = df[df["atg_to_length(%)"] >= perc]
    else:
        pcr = df[df["seq_length"] >= min_amino_length]

    raw_seq = pcr["sequence"]
    
    global_aa_counter = Counter()

    # 2. Iterate, Translate, and Profile
    for original_index, seq_value in raw_seq.items():
        if pd.isna(seq_value): continue
            
        if molecule_type == "DNA":
            seq_obj = Seq(seq_value[:len(seq_value) - len(seq_value)%3])
            chunks = str(seq_obj.translate(table=tb)).split("*")
        elif molecule_type == "RNA":
            seq_obj = Seq(seq_value[:len(seq_value) - len(seq_value)%3])
            chunks = str(seq_obj.back_transcribe().translate(table=tb)).split("*")
        else:
            chunks = [str(seq_value)]

        for chunk in chunks:
            # 1. Isolate the valid ORF within the chunk
            if molecule_type in ["DNA", "RNA"]:
                m_index = chunk.find("M")
                if m_index != -1:
                    valid_protein = chunk[m_index:]
                else:
                    valid_protein = "" 
            else:
                valid_protein = chunk

            amino_length = len(valid_protein)
            
            # 2. Enforce localized length and calculate biometrics
            if amino_length >= min_amino_length:
                gravy, mw, iso_p = None, None, None
                
                try:
                    clean_seq = valid_protein.replace("X", "").replace("B", "").replace("Z", "")
                    if len(clean_seq) > 0:
                        pa = ProteinAnalysis(clean_seq)
                        gravy = round(pa.gravy(), 3)
                        mw = round(pa.molecular_weight() / 1000, 2) 
                        iep = IsoelectricPoint(clean_seq)
                        iso_p = round(iep.pi(max_=14.0),2)
                except Exception:
                    pass 
                    
                extracted_proteins.append({
                    "original_index": original_index,
                    "protein_label": f"Protein_{count}",
                    "sequence": valid_protein,
                    "length": amino_length,
                    "MW (kDa)": mw,
                    "GRAVY Score": gravy,
                    "Isoelectric Point": iso_p
                })
                count += 1
                
                global_aa_counter.update(valid_protein)

    # 3. Compile the Protein DataFrame
    c_genes_df = pd.DataFrame(extracted_proteins)
    if not c_genes_df.empty:
        c_genes_df.set_index("original_index", inplace=True)

    # 4. Process the Global Counter into percentages
    total_aa = sum(global_aa_counter.values())
    if total_aa > 0:
        aa_data = [{"Amino Acid": aa, "Concentration (%)": round((cnt / total_aa) * 100, 2)} 
                   for aa, cnt in global_aa_counter.items()]
        final_composition_df = pd.DataFrame(aa_data).sort_values("Concentration (%)", ascending=False).reset_index(drop=True)
        final_composition_df.index += 1
    else:
        final_composition_df = pd.DataFrame(columns=["Amino Acid", "Concentration (%)"])
        
    return final_composition_df, c_genes_df