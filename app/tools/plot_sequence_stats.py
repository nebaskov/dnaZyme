import os

import pandas as pd
import plotly.express as px
from dotenv import load_dotenv

from tools import Counter

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')
PICTURES_PATH = os.getenv('PICTURES_PATH')


def get_sequences(path: str) -> list[str]:
    if path.endswith('xlsx'):
        df = pd.read_excel(path)
    elif path.endswith('csv'):
        df = pd.read_csv(path)
    else:
        raise ValueError('Not supported file format.')

    sequences: list[str] = df.sequence.tolist()
    return sequences


def count_tokens(
    sequences: list[str]
) -> tuple[dict[int, dict[str, int]], list[str]]:
    counter = Counter()
    count_dict: dict[int, dict[str, int]] = counter.count_in_sequences(
        sequences
    )
    return count_dict, counter.vocabulary


def count_dict_to_dataframe(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str]
) -> pd.DataFrame:
    count_df = pd.DataFrame(counts).transpose()
    count_df['position'] = count_df.index.tolist()

    for token in vocabulary:
        count_df[token] /= count_df[token].sum()

    return count_df


def plot_histograms(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str]
) -> None:
    count_df: pd.DataFrame = count_dict_to_dataframe(counts, vocabulary)
    fig = px.bar(
        count_df,
        x='position',
        y=vocabulary,
        title='Nucleotide position-wise frequency',
        labels={
            'variable': 'Nucleotide',
            'value': 'frequency'
        },
    )
    fig.show()
    if input('Confirm image saving: y/n? ') == 'y':
        fig.write_image(
            os.path.join(PICTURES_PATH, 'sequences_barplot.png'),
            scale=4
        )
    else:
        print('No image saving')


def plot_lines(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str]
) -> None:
    count_df: pd.DataFrame = count_dict_to_dataframe(counts, vocabulary)
    fig = px.line(
        count_df,
        x='position',
        y=vocabulary,
        title='Nucleotide position-wise frequency',
        labels={
            'variable': 'Nucleotide',
            'value': 'frequency'
        },
    )
    fig.show()
    if input('Confirm image saving: y/n? ') == 'y':
        fig.write_image(
            os.path.join(PICTURES_PATH, 'sequences_lineplot.png'),
            scale=4
        )
    else:
        print('No image saving')


def main(path: str) -> None:
    sequences: list[str] = get_sequences(path)
    count_dict, vocabulary = count_tokens(sequences)
    # plot_histograms(count_dict, vocabulary)
    plot_lines(count_dict, vocabulary)


if __name__ == '__main__':
    main(path=os.path.join(DATA_PATH, 'aligned_sequences.csv'))
