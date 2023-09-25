import os
from dotenv import load_dotenv

import pandas as pd
import numpy as np
from sklearn.manifold import MDS
from Levenshtein import distance
from sklearn.neighbors import kneighbors_graph
from sklearn.cluster import AgglomerativeClustering

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from tools import count_tokens, count_dict_to_dataframe

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')
PICTURES_PATH = os.getenv('PICTURES_PATH')


def get_levenshtein_analysis(sequence_df: pd.DataFrame) -> tuple[
    AgglomerativeClustering,
    np.array
]:
    all_sequences = list(dict.fromkeys(sequence_df.e.tolist()))

    distances = np.zeros((len(all_sequences), len(all_sequences)))
    for i in range(len(all_sequences)):
        for j in range(i + 1, len(all_sequences)):
            distances[i, j] = distance(all_sequences[i], all_sequences[j])
            distances[j, i] = distances[i, j]

    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=1)
    pos = mds.fit_transform(distances)

    k = 2
    connectivity = kneighbors_graph(pos, n_neighbors=k, include_self=False)
    model = AgglomerativeClustering(
        n_clusters=k,
        linkage='ward',
        connectivity=connectivity
    )
    model.fit(pos)
    return model, pos


def plot_levenshtein_analysis_1(sequence_df: pd.DataFrame):
    sequence_df.drop_duplicates(subset='e', inplace=True)
    model, pos = get_levenshtein_analysis(sequence_df)
    log_kobs = np.log10(sequence_df['kobs'])
    norm = plt.Normalize(log_kobs.min(), log_kobs.max())
    cmap = ListedColormap([
        "#2E4451",
        "#9C5551",
        "#747178",
        "#E0B3AD",
        "#DAD1D2"
    ])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    ax1 = sns.scatterplot(
        x=pos[:, 0],
        y=pos[:, 1],
        c=log_kobs,
        alpha=0.7,
        cmap=cmap
    )
    ax1.figure.colorbar(sm, ax=ax1, label='log(kobs)')
    return ax1


def plot_levenshtein_analysis_2(sequence_df: pd.DataFrame):
    sequence_df.drop_duplicates(subset='e', inplace=True)
    model, pos = get_levenshtein_analysis(sequence_df)
    cmap = ListedColormap(["#2E4451", "#9C5551"])

    ax2 = sns.scatterplot(
        x=pos[:, 0],
        y=pos[:, 1],
        c=model.labels_,
        cmap=cmap
    )

    return ax2


def plot_frequency(sequence_df: pd.DataFrame):
    count_dict, vocabulary = count_tokens(sequence_df.sequence.tolist())
    count_df = count_dict_to_dataframe(count_dict, vocabulary, scaled=True)
    token_df = count_df[vocabulary]
    # colors: list[str] = [
    #     '#cc00ff',  # purple
    #     '#ff5050',  # green
    #     '#00cc99',  # red
    #     '#3366ff',  # blue
    #     # '#ffff99'  # yellow
    # ]
    colors = [
        "#2E4451",
        "#9C5551",
        "#747178",
        "#E0B3AD"
    ]

    return token_df, colors


def plot_full_analysis(
        sequence_df: pd.DataFrame,
        aligned_seq: pd.DataFrame
):
    fig, ax = plt.subplots(2, 2, dpi=200)
    plt.sca(ax[1, 0])
    plot_levenshtein_analysis_1(sequence_df)

    plt.sca(ax[0, 0])
    ax[0, 0].set_title('DNAzyme sequence clustering')
    plot_levenshtein_analysis_2(sequence_df)

    ax[0, 1].set_title('DNAzyme sequence alignment')

    plt.sca(ax[1, 1])
    token_df, colors = plot_frequency(aligned_seq)
    token_df.plot(
        kind='bar',
        stacked=True,
        color=colors,
        ylabel='Percentage',
        ax=ax[1, 1],
        rot='horizontal'
    )
    ax[1, 1].legend(bbox_to_anchor=(1.1, 0.8), loc='right')
    plt.setp(ax[1, 1].get_xticklabels(), visible=False)

    fig.set_figheight(14)
    fig.set_figwidth(25)
    fig.tight_layout()
    
    # set the spacing between subplots
    plt.subplots_adjust(
        left=0.1,
        bottom=0.1,
        right=0.9,
        top=0.9,
        wspace=0.4,
        hspace=0.4
    )

    fig.savefig(
        os.path.join(PICTURES_PATH, 'sequence_analysis.png')
    )


if __name__ == '__main__':
    seq_df = pd.read_excel(
        os.path.join(DATA_PATH, 'preprocessed_dataset.xlsx')
    )
    align_df = pd.read_csv(
        os.path.join(DATA_PATH, 'clean_align_trend_seq.csv')
    )
    plot_full_analysis(seq_df, align_df)
