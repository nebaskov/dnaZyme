import os

import pandas as pd
from Bio import SeqIO
from dotenv import load_dotenv

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')


def read_fasta(path: str) -> dict[str, str]:
    fasta_sequences = SeqIO.parse(open(path), 'fasta')
    dict_sequences: dict[str, str] = {
        fasta.id: str(fasta.seq) for fasta in fasta_sequences
    }
    return dict_sequences


def fasta_to_dataframe(
    read_path: str,
    save_path: str = None,
    save: bool = False
) -> pd.DataFrame:
    processed_sequences: dict[str, str] = read_fasta(read_path)

    data = pd.DataFrame(columns=['id', 'sequence'])
    data['id'] = list(processed_sequences.keys())
    data['sequence'] = list(processed_sequences.values())

    if save and save_path is not None:
        data.to_csv(save_path, index=False)

    return data


if __name__ == '__main__':
    filepath = 'aligned_sequences.fasta'
    processed_sequences: dict[str, str] = read_fasta(
        os.path.join(DATA_PATH, filepath)
    )
    data = pd.DataFrame(columns=['id', 'sequence'])
    data['id'] = list(processed_sequences.keys())
    data['sequence'] = list(processed_sequences.values())
    data.to_csv(
        os.path.join(DATA_PATH, 'aligned_sequences.csv'),
        index=False
    )
