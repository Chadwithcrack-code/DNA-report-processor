from Bio import SeqIO
import pandas as pd
import io



def process(open_file):
    try:
        # Wrap the byte stream in a TextIOWrapper for memory-efficient text streaming.
        # This prevents loading the entire 50MB+ file into memory at once.
        text_stream = io.TextIOWrapper(open_file, encoding="utf-8")
        
        # This yields one record at a time directly into the DataFrame engine.
        record_generator = (
            {
                "id": record.id,
                "description": record.description,
                "sequence": str(record.seq)
            }
            for record in SeqIO.parse(text_stream, "fasta")
        )
        
        # Consume the generator lazily
        df = pd.DataFrame.from_records(record_generator)
        print(df)
        
        if df.empty:
            raise ValueError("The uploaded file contains no valid FASTA sequences.")
            
        return df
        
    except UnicodeDecodeError:
        raise ValueError("Encoding error: Ensure the uploaded file is an unzipped, valid UTF-8 FASTA file.")
    except Exception as e:
        # Catch any unexpected BioPython parsing errors
        raise RuntimeError(f"Failed to parse FASTA file: {str(e)}")
    finally:
        # Reset the file pointer so Streamlit can read it again if the user changes a parameter
        open_file.seek(0)
