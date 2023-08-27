import os

import pandas as pd
from dotenv import load_dotenv

from tools import convert_to_fasta

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')

data = pd.read_excel(DATA_PATH + 'preprocessed_dataset.xlsx')
convert_to_fasta(data, 'e', 'id', DATA_PATH + 'sequences.fasta')
