import pandas as pd


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
