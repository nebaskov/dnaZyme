import os
from typing import Literal

import pandas as pd
from dotenv import load_dotenv

from scipy.stats import ttest_1samp, ttest_rel, ttest_ind, chi2_contingency, chi2

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from tools import Counter

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')
PICTURES_PATH = os.getenv('PICTURES_PATH')


def index_generator(data: pd.Index, chunk_size: int) -> pd.Index:
    for i in range(0, len(data), chunk_size):
        yield data[i:i+chunk_size]


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
    for position, tokens in count_dict.items():
        tokens['token_sum'] = sum(tokens.values())
        count_dict[position] = tokens
    return count_dict, counter.vocabulary


def count_dict_to_dataframe(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str],
    scaled: bool = False
) -> pd.DataFrame:
    count_df = pd.DataFrame(counts).transpose()
    count_df['position'] = count_df.index.tolist()

    if scaled:
        for token in vocabulary:
            count_df[token] /= count_df['token_sum']
            count_df[token] *= 100

    return count_df


def get_sample_distribution(
    count_df: pd.DataFrame,
    vocabulary: list[str]
) -> pd.DataFrame:
    sample_df = count_df[['position', 'token_sum']].copy()
    sample_df['sample_value'] = sample_df['token_sum'].apply(
        lambda x: int(round(x * (1 / len(vocabulary))))
    )
    return sample_df


def run_ttest(
    count_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    vocabulary: list[str]
) -> pd.DataFrame:
    p_values: list[float] = []
    count_df_copy = count_df.copy()
    df = 55
    for idx in count_df.index:
        cdf = chi2.cdf(count_df.iloc[idx][vocabulary].tolist(), df)
        p_val = 2 * (1 - cdf)
        chi2, p_val, dof, expected = chi2_contingency(
            [
                count_df.iloc[idx][vocabulary].tolist(),
                [0.25, 0.25, 0.25, 0.25]
            ]
            # sample_df['sample_value'][idx]
        )
        p_values.append(p_val)
    count_df_copy['p_value'] = p_values
    return count_df_copy


def plot_ttest_histogram(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str],
) -> None:
    count_df: pd.DataFrame = count_dict_to_dataframe(
        counts=counts,
        vocabulary=vocabulary,
        scaled=True
    )
    count_df['position'] += 1
    sample_df: pd.DataFrame = get_sample_distribution(count_df, vocabulary)
    ttest_df: pd.DataFrame = run_ttest(count_df, sample_df, vocabulary)

    subplots = make_subplots(
        rows=2,
        cols=2,
        shared_yaxes='all',
        row_width=[0.4, 0.4],
        vertical_spacing=0.1,
        horizontal_spacing=0.05
    )

    chunk_size = 34
    n_chunks: int = ttest_df.shape[0] // chunk_size + \
        (1 if ttest_df.shape[0] % chunk_size != 0 else 0)

    colors: list[str] = [
        '#cc00ff',  # purple
        '#ff5050',  # green
        '#00cc99',  # red
        '#3366ff',  # blue
        '#ffff99'   # yellow
    ]
    colors_dict: dict[str, str] = dict(
        zip(vocabulary, colors[:len(vocabulary)])
    )

    for token in sorted(vocabulary, reverse=False):
        data_iterator = index_generator(ttest_df.index, chunk_size)
        for n in range(n_chunks):
            idx = next(data_iterator)
            subplots.add_trace(
                go.Bar(
                    name=token,
                    x=ttest_df['position'][idx],
                    y=ttest_df.loc[idx, token],
                    legendgroup=token,
                    marker_color=colors_dict.get(token),
                    showlegend=(True if n == 0 else False),
                ),
                (n // 2 + 1),
                (n % 2 + 1)
            )
    subplots.update_layout(
        title='Nucleotide position-wise frequency',
        barmode='stack',
        height=800,
        width=1200
    )
    subplots.update_xaxes(title='Position')
    subplots['layout']['yaxis1']['title'] = 'Frequency'
    subplots['layout']['yaxis3']['title'] = 'Frequency'
    subplots.show()
    if input('Confirm image saving: y/n? ') == 'y':
        subplots.write_image(
            os.path.join(PICTURES_PATH, 'sequences_barplot_no_norm.png'),
            scale=4
        )
    else:
        print('No image saving')


def plot_histograms(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str],
    scaled: bool = False
) -> None:
    count_df: pd.DataFrame = count_dict_to_dataframe(
        counts=counts,
        vocabulary=vocabulary,
        scaled=scaled
    )
    count_df['position'] += 1

    subplots = make_subplots(
        rows=2,
        cols=2,
        shared_yaxes='all',
        row_width=[0.4, 0.4],
        vertical_spacing=0.1,
        horizontal_spacing=0.05
    )

    chunk_size = 50
    n_chunks: int = count_df.shape[0] // chunk_size + \
        (1 if count_df.shape[0] % chunk_size != 0 else 0)

    colors: list[str] = [
        '#cc00ff',  # purple
        '#ff5050',  # green
        '#00cc99',  # red
        '#3366ff',  # blue
        '#ffff99'   # yellow
    ]
    colors_dict: dict[str, str] = dict(
        zip(vocabulary, colors[:len(vocabulary)])
    )

    for token in sorted(vocabulary, reverse=False):
        data_iterator = index_generator(count_df.index, chunk_size)
        for n in range(n_chunks):
            idx = next(data_iterator)
            subplots.add_trace(
                go.Bar(
                    name=token,
                    x=count_df['position'][idx],
                    y=count_df.loc[idx, token],
                    legendgroup=token,
                    marker_color=colors_dict.get(token),
                    showlegend=(True if n == 0 else False),
                ),
                (n // 2 + 1),
                (n % 2 + 1)
            )
    subplots.update_layout(
        title='Nucleotide position-wise frequency',
        barmode='stack',
        height=800,
        width=1200
    )
    subplots.update_xaxes(title='Position')
    subplots['layout']['yaxis1']['title'] = 'Frequency'
    subplots['layout']['yaxis3']['title'] = 'Frequency'
    subplots.show()
    if input('Confirm image saving: y/n? ') == 'y':
        subplots.write_image(
            os.path.join(PICTURES_PATH, 'sequences_barplot_no_norm.png'),
            scale=4
        )
    else:
        print('No image saving')


def plot_lines(
    counts: dict[int, dict[str, int]],
    vocabulary: list[str],
    scaled: bool = False
) -> None:
    count_df: pd.DataFrame = count_dict_to_dataframe(
        counts=counts,
        vocabulary=vocabulary,
        scaled=scaled
    )
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


def plot_seq_stats(
    filepath: str,
    options: list[
        Literal['barplot'] | Literal['lineplot'] | Literal['ttest_histogram']
    ]
) -> None:
    sequences: list[str] = get_sequences(filepath)
    count_dict, vocabulary = count_tokens(sequences)
    if 'barplot' in options:
        plot_histograms(count_dict, vocabulary)
    if 'lineplot' in options:
        plot_lines(count_dict, vocabulary)
    if 'ttest_histogram' in options:
        plot_ttest_histogram(count_dict, vocabulary)



if __name__ == '__main__':
    plot_seq_stats(
        filepath=os.path.join(DATA_PATH, 'clean_align_trend_seq.csv'),
        options=[
            'barplot',
            # 'lineplot',
            # 'ttest_histogram'
        ]
    )
