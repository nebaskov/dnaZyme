import os

import pandas as pd
import plotly.express as px
from dotenv import load_dotenv

from tools import Counter

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')


def get_sequences(path: str) -> list[str]:
    df = pd.read_excel(path)
    sequences: list[str] = df.e.tolist()
    return sequences


def count_tokens(
    sequences: list[str]
) -> tuple[dict[int, dict[str, int]], list[str]]:
    counter = Counter()
    count_dict: dict[int, dict[str, int]] = counter.count_in_sequences(
        sequences
    )
    return count_dict, counter.vocabulary


def plot_histograms(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str]
) -> None:
    count_df = pd.DataFrame(counts).transpose()
    count_df['position'] = count_df.index.tolist()
    fig = px.bar(
        count_df,
        x='position',
        y=vocabulary,
        title='Nucleotide position-wise frequency',
        labels={
            'variable': 'Nucleotide',
            'value': 'frequency'
        }
    )
    fig.show()
    if input('Confirm image saving: y/n') == 'y':
        fig.write_image('pictures/sequences_barplot.png')
    else:
        print('No image saving')


def main(path: str) -> None:
    sequences: list[str] = get_sequences(path)
    count_dict, vocabulary = count_tokens(sequences)
    plot_histograms(count_dict, vocabulary)


if __name__ == '__main__':
    main(path=os.path.join(DATA_PATH, 'preprocessed_dataset.xlsx'))
