import pandas as pd
from Bio import SeqIO


def convert_to_fasta(
    seq_df: pd.DataFrame,
    seq_column: str | pd.Index,
    id_column: str | pd.Index,
    path: str,
) -> None:
    sequences: list[str] = seq_df[seq_column].tolist()
    ids: list[int] = seq_df[id_column].tolist()
    fasta = open(path, 'w')
    for idx, seq in zip(ids, sequences):
        fasta.write('>' + str(idx) + '\n' + seq + '\n')
    fasta.close()


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


def process_fasta(read_path: str, save_path: str) -> None:
    raw_fasta: pd.DataFrame = fasta_to_dataframe(read_path)
    raw_fasta.drop_duplicates(subset='sequence', inplace=True)
    raw_fasta.sort_values(by='id', inplace=True)
    convert_to_fasta(
        raw_fasta,
        'sequence',
        'id',
        save_path
    )