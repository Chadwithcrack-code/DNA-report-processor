from Bio import SeqIO
from io import StringIO
import pandas as pd
import re

def process(open_file):
    file = StringIO(open_file.getvalue().decode("utf-8"))
    data = []

    print(file)
    print("this line of code is running")

    for record in SeqIO.parse(file,"fasta"):
            data.append({
                "id": record.id,
                "description": record.description,
                "sequence": str(record.seq)
            })
            print(record)

    return pd.DataFrame.from_records(data)

