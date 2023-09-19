import os
from dotenv import load_dotenv

from tools import process_fasta, fasta_to_dataframe

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')


if __name__ == '__main__':
    # process_fasta(
    #     read_path=os.path.join(DATA_PATH, 'sequences.fasta'),
    #     save_path=os.path.join(DATA_PATH, 'clean_seq.fasta')
    # )
    fasta_to_dataframe(
        read_path=os.path.join(DATA_PATH, 'clean_align_trend_seq.fasta'),
        save_path=os.path.join(DATA_PATH, 'clean_align_trend_seq.csv'),
        save=True
    )