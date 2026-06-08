import streamlit as st
from processor import process
from cleaner import clean
import analysis

if "counter" not in st.session_state:
    st.session_state.counter = 0
    st.session_state.files = ""


def main():
    st.title("CO Automation FASTA PROCESSOR")
    if st.session_state.counter == 0:
        st.header("Upload file below")
        st.session_state.files = st.file_uploader("Insert",type=[".fasta",".fa",".fas",".fna"])
    
    if st.session_state.files != None:
        df = clean(process(st.session_state.files))
        amino = analysis.analyze_amino(df=df)
        st.write(df)
        st.write(amino)
        st.write(df.columns)

    


if __name__ == "__main__":
    main()
