import streamlit as st
from processor import process
from cleaner import clean
import io
import analysis
import traceback
import matplotlib.pyplot as plt


st.set_page_config(page_title="FASTA Processor", layout="wide")

# Cache
@st.cache_data(show_spinner=False, max_entries=5, ttl=3600)
def get_clean_df(file_bytes):
    file_obj = io.BytesIO(file_bytes)
    v = clean(process(file_obj))
    return v



def main():
    st.title("🧬 CO Automation FASTA PROCESSOR")
    st.markdown(
        """
        <div style='margin-top: -15px; margin-bottom: 20px; color: gray;'>
            <em>Created by <a href='https://www.linkedin.com/in/chad-o-ntyamba-487779413/' target='_blank' style='color: #4CAF50; text-decoration: none;'>Chad Ontyamba</a> | CO Automation</em>
        </div>
        """, 
        unsafe_allow_html=True
    )
    
    st.markdown("Automated sequence profiling, translation, and biophysical analysis.")
    st.divider()


    # SIDEBAR PARAMETERS

    st.sidebar.title("Parameters")
    st.sidebar.divider()
    
    min_length = st.sidebar.number_input("Min protein length", value=50, step=10)
    st.sidebar.divider()
    
    threshold = st.sidebar.slider("Start codon threshold(%)", min_value=0.01, max_value=10.0, value=3.0, step=0.1)
    st.sidebar.divider()
    
    Translation_t = {
        "Standard (Table 1)": 1,
        "Bacterial/Archaeal (Table 11)": 11,
        "Mitochondrial (Table 2)": 2,
        "Yeast Mitochondrial (Table 3)": 3
    }
    # Extract the integer ID for the backend logic
    selected_table = st.sidebar.selectbox("Translation table", options=list(Translation_t.keys()))
    table_id = Translation_t[selected_table]

    # FILE INGESTION:

    st.header("1. Upload Data")
    uploaded_file = st.file_uploader("Upload FASTA (Max 200MB)", type=[".fasta", ".fa", ".fas", ".fna"])

    # Streamlit naturally waits here until a file is provided. No counter needed.
    if uploaded_file is not None:
        
        st.header("2. Analysis Report")
        with st.spinner("Processing dataset..."):
            
            try:
                # 1. Safely extract bytes and pass to cached ingestion
                file_bytes = uploaded_file.getvalue()
                df = get_clean_df(file_bytes)
                
                # 2. Run your analysis algorithm
                amino, proteins = analysis.analyze_amino(
                    df=df,
                    perc=threshold,
                    min_amino_length=min_length,
                    tb=table_id
                )
                
                # 3. Dashboard UI
                st.success("✅ Analysis complete!")
                
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("Total Sequences Parsed", value=len(df))

                if not amino.empty:
                    top_index = amino["Concentration (%)"].idxmax()
                    top_aa = amino.loc[top_index, "Amino Acid"]
                    top_conc = amino.loc[top_index, "Concentration (%)"]
                    col2.metric(label=f"Most Expressed: {top_aa}", value=f"{top_conc:.2f}%")
                else:
                    col2.metric(label="Most Expressed", value="None found")
                    
                col3.metric("Longest Polypeptide", value=proteins['sequence'].str.len().max() if not proteins.empty else 0)
                
                with col4:
                    fig4, ax4 = plt.subplots(figsize=(10, 6))
                    ax4.hist(df["GC (%)"], bins=20, color='skyblue',edgecolor='black')
                    ax4.set_xlabel("GC (%)")
                    ax4.set_ylabel("Number of Sequences")
                    ax4.set_title("GC (%) Distribution")
                    st.pyplot(fig4)
                    plt.close(fig4)
                    


                # --- NEW METRICS BLOCK ---
                if not proteins.empty and "MW (kDa)" in proteins.columns:
                    st.divider()
                    col4, col5, col6 = st.columns(3)
                    avg_mw = proteins["MW (kDa)"].mean()
                    avg_gravy = proteins["GRAVY Score"].mean()
                    avg_iso = proteins["Isoelectric Point"].mean()
                    
                    col4.metric("Avg Molecular Weight", f"{avg_mw:.2f} kDa", help="Calculated using BioPython ProtParam")
                    col5.metric("Avg GRAVY Score", f"{avg_gravy:.2f}", help="Grand Average of Hydropathy. Positive values indicate hydrophobic proteins.")
                    col6.metric("Avg Isoelectric Point", f"{avg_iso:.2f}", help="The pH at which the molecule carries no net electrical charge.")
                # -------------------------

                st.divider()
    
                st.subheader("Protein Profile")
                
                # Only plot the chart if there is data
                if not amino.empty:
                    col7, col8, col9 = st.columns(3)
                    
                    with col7:

                        fig1, ax = plt.subplots(figsize=(10, 6))
                        ax.scatter(proteins["Isoelectric Point"], proteins["MW (kDa)"])
                        ax.set_xlabel("Isoelectric Point")
                        ax.set_ylabel("Molecular Weight (kDa)")
                        ax.set_title("Virtual 2D Gel (pI vs. Molecular Weight)")
                        st.pyplot(fig1)
                        plt.close(fig1)

                    with col8:
                        fig2, ax2 = plt.subplots(figsize=(10, 6))
                        ax2.hexbin(proteins["Isoelectric Point"], proteins["GRAVY Score"], gridsize=50, cmap='inferno')
                        ax2.set_xlabel("Isoelectric Point")
                        ax2.set_ylabel("GRAVY Score")
                        ax2.set_title("Proteome Hydrophobicity Landscape (GRAVY vs. pI)")
                        st.pyplot(fig2)
                        plt.close(fig2)

                    with col9:
                        fig3, ax3 = plt.subplots(figsize=(10, 6))
                        ax3.bar(amino["Amino Acid"], amino["Concentration (%)"])
                        ax3.set_xlabel("Amino Acid")
                        ax3.set_ylabel("Concentration (%)")
                        ax3.set_title("Amino Acid Concentration")
                        st.pyplot(fig3)
                        plt.close(fig3)
                else:
                    st.info("No amino acids passed the filtering thresholds. Try lowering the Start Codon Threshold or Minimum Protein Length.")
                
                st.subheader("Extracted Coding Sequences (CDS) - Preview")
                
                #  UI SAFETY LIMIT 
                # Only render the first 1000 rows to the browser to prevent UI crashing.
                # The browser stays light and fast, no matter how big the genome is.
                max_display_rows = 500
                st.dataframe(proteins.head(max_display_rows), use_container_width=True)
                
                if len(proteins) > max_display_rows:
                    st.info(f"UI Limited to the first {max_display_rows} sequences to maintain browser performance. The full dataset contains {len(proteins):,} sequences.")

                # The download button still exports the FULL dataset, not just the preview.
                st.download_button(
                    label="📥 Export Full Protein CDS (CSV)",
                    data=proteins.to_csv(index=False).encode('utf-8'), 
                    file_name="co_automation_proteins_report.csv",
                    mime="text/csv",
                    use_container_width=True
                )

            except Exception as e:
                # Global exception trap catches processor, cleaner, and analysis errors
                st.error(f"Pipeline Interrupted: {str(e)}")
                # Expandable traceback is excellent for debugging but keeps the UI clean for the client
                with st.expander("Show System Logs"):
                    st.code(traceback.format_exc(), language="bash")

if __name__ == "__main__":
    main()