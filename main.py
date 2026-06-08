import streamlit as st
from processor import process
from cleaner import clean
from io import BytesIO
import analysis

if "counter" not in st.session_state:
    st.session_state.counter = 0
    st.session_state.files = ""
    st.session_state.min_amino_length = 50
    st.session_state.start_codon_threshold = 2.8
    st.session_state.clean_df = ""
    st.session_state.translation = 1


def control(step):
    st.session_state.counter = step

@st.cache_data(show_spinner=False)
def get_clean_df(file_obj):
    v = clean(process(file_obj))
    return v


def main():
    st.set_page_config(page_title="FASTA Processor", layout="wide")
    st.sidebar.title("Parameters")
    st.sidebar.divider()
    st.sidebar.number_input("Min protein length", value=50, key="min_amino_length")
    st.sidebar.divider()
    st.sidebar.slider("Start codon threshold(%)", value=3, key="start_codon_threshold")
    st.sidebar.divider()
    
    Translation_t = {
    "Standard (Table 1)": 1,
    "Bacterial/Archaeal (Table 11)": 11,
    "Mitochondrial (Table 2)": 2,
    "Yeast Mitochondrial (Table 3)": 3
    }
    st.session_state.translation = st.sidebar.selectbox("Translation table", options=list(Translation_t.keys()))


    if st.session_state.counter == 0:
        st.title("🧬 CO Automation FASTA PROCESSOR")
        st.markdown("Automated sequence profiling, translation, and biophysical analysis.")


        st.header("Upload file below")
        st.session_state.files = st.file_uploader("Insert",type=[".fasta",".fa",".fas",".fna"])
        st.button("Process", on_click=control, args=(1,))


    with st.spinner("Processing data... Please wait."):
    
        if st.session_state.files is not None and st.session_state.counter == 1:
            df = get_clean_df(st.session_state.files)

            
            try:
                amino, proteins = analysis.analyze_amino(df=df,
                perc=st.session_state.start_codon_threshold,
                min_amino_length=st.session_state.min_amino_length,
                tb=st.session_state.translation

                )
                st.success("Analysis complete!!")
                st.metric("Longest poly-peptide chain:", value=proteins['sequence'].str.len().max())
                st.success("Analysis successful!")
                col1, col2, col3 = st.columns(3)
                col1.metric("Total Sequences Parsed", value=len(df))
                #st.write(amino.loc[amino["Concentration (%)"].idxmax()])
                # 1. Find the index of the highest concentration
                top_index = amino["Concentration (%)"].idxmax()

                # 2. Extract the pure string and float values using that index
                top_aa = amino.loc[top_index, "Amino Acid"]
                top_conc = amino.loc[top_index, "Concentration (%)"]

                # 3. Pass the clean variables into the Streamlit metric
                col2.metric(label=f"Most Expressed Amino acid: {top_aa}", value=f"{top_conc:.2f}%")
                col3.metric("Longest Polypeptide Chain", value=proteins['sequence'].str.len().max() if not proteins.empty else 0)

                st.divider()
    
                st.subheader("Amino Acid Distribution Profile")
                st.bar_chart(data=amino, x="Amino Acid", y="Concentration (%)", color="#4CAF50", horizontal=True)
                st.subheader("Extracted Coding Sequences (CDS)")
                st.dataframe(proteins, use_container_width=True)

                st.download_button(
                    label="📥 Export Protein CDS (CSV)",
                    data=proteins.to_csv().encode('utf-8'),
                    file_name="cleaned_proteins_report.csv",
                    mime="text/csv",
                    use_container_width=True
                )

            except ValueError as e:
                st.error(str(e))
                st.write(df)
    

    


if __name__ == "__main__":
    main()
