import streamlit as st
from processor import process
from cleaner import clean
import io
import analysis
import traceback

# Setup must be the first Streamlit command
st.set_page_config(page_title="FASTA Processor", layout="wide")

# Cache using the raw bytes to prevent Streamlit UploadedFile object dropping
@st.cache_data(show_spinner=False)
def get_clean_df(file_bytes):
    file_obj = io.BytesIO(file_bytes)
    v = clean(process(file_obj))
    return v

def main():
    st.title("🧬 CO Automation FASTA PROCESSOR")
    st.markdown("Automated sequence profiling, translation, and biophysical analysis.")

    # ==========================================
    # SIDEBAR PARAMETERS
    # ==========================================
    st.sidebar.title("Parameters")
    st.sidebar.divider()
    
    min_length = st.sidebar.number_input("Min protein length", value=50, step=10)
    st.sidebar.divider()
    
    threshold = st.sidebar.slider("Start codon threshold(%)", min_value=0.5, max_value=10.0, value=3.0, step=0.1)
    st.sidebar.divider()
    
    Translation_t = {
        "Standard (Table 1)": 1,
        "Bacterial/Archaeal (Table 11)": 11,
        "Mitochondrial (Table 2)": 2,
        "Yeast Mitochondrial (Table 3)": 3
    }
    # Extract the actual integer ID for the backend logic
    selected_table = st.sidebar.selectbox("Translation table", options=list(Translation_t.keys()))
    table_id = Translation_t[selected_table]

    # ==========================================
    # FILE INGESTION
    # ==========================================
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
                
                col1, col2, col3 = st.columns(3)
                col1.metric("Total Sequences Parsed", value=len(df))
                
                top_index = amino["Concentration (%)"].idxmax()
                top_aa = amino.loc[top_index, "Amino Acid"]
                top_conc = amino.loc[top_index, "Concentration (%)"]
                
                col2.metric(label=f"Most Expressed: {top_aa}", value=f"{top_conc:.2f}%")
                col3.metric("Longest Polypeptide", value=proteins['sequence'].str.len().max() if not proteins.empty else 0)

                st.divider()
    
                st.subheader("Amino Acid Distribution Profile")
                st.bar_chart(data=amino, x="Amino Acid", y="Concentration (%)", color="#4CAF50", horizontal=True)
                
                st.subheader("Extracted Coding Sequences (CDS)")
                st.dataframe(proteins, use_container_width=True)

                st.download_button(
                    label="📥 Export Protein CDS (CSV)",
                    data=proteins.to_csv(index=False).encode('utf-8'), # Removed index for cleaner client exports
                    file_name="cleaned_proteins_report.csv",
                    mime="text/csv",
                    use_container_width=True
                )

            except Exception as e:
                # Global exception trap catches processor, cleaner, and analysis errors
                st.error(f"⚠️ Pipeline Interrupted: {str(e)}")
                # Expandable traceback is excellent for debugging but keeps the UI clean for the client
                with st.expander("Show System Logs"):
                    st.code(traceback.format_exc(), language="bash")

if __name__ == "__main__":
    main()