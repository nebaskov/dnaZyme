import os
from dotenv import load_dotenv
from biotite.sequence import graphics, NucleotideSequence
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from tools import read_fasta

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')
PICTURES_PATH = os.getenv('PICTURES_PATH')


def plot_similarity(path: str) -> None:
    seq: list[str] = list(read_fasta(path).values())
    seq_obj = [NucleotideSequence(s, ambiguous=True) for s in seq]
    _, ax = plt.subplots(1, 1)
    graphics.plot_alignment_similarity_based(axes=ax, alignment=seq_obj)


if __name__ == '__main__':
    path = os.path.join(DATA_PATH, 'aligned_sequences.fasta')
    plot_similarity(path)
