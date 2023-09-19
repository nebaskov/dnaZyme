import os

import pandas as pd
from dotenv import load_dotenv
import biotite.sequence as seq
from biotite.application.clustalo.app import ClustalOmegaApp, Alignment

from tools import read_fasta

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')


def read_file(path: str, seq_column: str = None):
    if path.endswith('.fasta'):
        sequences: list[str] = list(read_fasta(path).values())
    elif path.endswith('.xlsx') and seq_column is not None:
        seq_df = pd.read_excel(path)
        sequences: list[str] = seq_df[seq_column].unique().tolist()
    elif path.endswith('.csv') and seq_column is not None:
        seq_df = pd.read_csv(path)
        sequences: list[str] = seq_df[seq_column].unique().tolist()
    else:
        raise ValueError('Not supported file format')
    return sequences


def align_sequences(path: str, seq_column: str = None) -> Alignment:
    """Perform multiple sequence alignment using Clustal-Omega tool.

    Args:
        path: path to fasta, excel or csv file containing sequences
        seq_column: the name of column containing sequence data if file type is excel or csv

    Returns:
        list[str]: list of aligned sequences
    """
    sequences: list[str] = read_file(path, seq_column)
    seq_classes: list[seq.NucleotideSequence] = [
        seq.NucleotideSequence(sequence) for sequence in sequences
    ]
    alignment = ClustalOmegaApp(seq_classes)
    return alignment


if __name__ == '__main__':
    path = os.path.join(DATA_PATH, 'preprocessed_dataset.xlsx')
    align = align_sequences(path, 'e')
