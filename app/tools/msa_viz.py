import os

from dotenv import load_dotenv
from pymsaviz import MsaViz

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')
PICTURES_PATH = os.getenv('PICTURES_PATH')

def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def sort_hamming_distance(s1):
    """Return the Hamming distance between equal-length sequences"""
    s2 = 'XCTGCAAGGAGAGGGTGCAGCGGGAGTGGGGGGTTGXTGGGTTCCAA'
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def visualize_msa(path: str) -> None:
    msa = MsaViz(
        msa=path,
        format='fasta',
        show_count=False,
        show_consensus=True
    )
    # msa.savefig(savefile=os.path.join(PICTURES_PATH, 'msa_viz.png'), dpi=200)
    consensus = msa.consensus_seq

    filtered = []
    for record in msa.seq_list:
        if hamming_distance(record, consensus) < 50:
            filtered.append(record)
    sorted_filt = sorted(filtered, key=sort_hamming_distance, reverse=False)

    with open('filted_seq.txt', 'w') as file:
        for seq in sorted_filt:
            file.write(seq + '\n')


if __name__ == '__main__':
    visualize_msa(os.path.join(DATA_PATH, 'clean_align_trend_seq.fasta'))
