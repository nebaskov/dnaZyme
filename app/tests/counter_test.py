from tools import Counter


def one_sequence_test() -> None:
    counter = Counter()
    sequences = ['AAT, GC']
    counter.count_in_sequences(sequences)
    counts_to_test = counter.count_dict
    sample_counts = {
        0: {'A': 1, 'T': 0, 'G': 0, 'C': 0},
        1: {'A': 1, 'T': 0, 'G': 0, 'C': 0},
        2: {'A': 0, 'T': 1, 'G': 0, 'C': 0},
        3: {'A': 0, 'T': 0, 'G': 1, 'C': 0},
        4: {'A': 0, 'T': 0, 'G': 0, 'C': 1},
    }
    assert counts_to_test == sample_counts, 'нихера не работает, переделывай'


def multiple_sequence_test() -> None:
    counter = Counter()
    sequences = ['AATGC', 'ATTG', 'GTTAAA']
    counter.count_in_sequences(sequences)
    counts_to_test = counter.count_dict
    sample_counts = {
        0: {'A': 2, 'T': 0, 'G': 1, 'C': 0},
        1: {'A': 1, 'T': 2, 'G': 0, 'C': 0},
        2: {'A': 0, 'T': 3, 'G': 0, 'C': 0},
        3: {'A': 1, 'T': 0, 'G': 2, 'C': 0},
        4: {'A': 1, 'T': 0, 'G': 0, 'C': 1},
        5: {'A': 1, 'T': 0, 'G': 0, 'C': 0},
    }
    assert counts_to_test == sample_counts, 'нихера не работает, переделывай'


def different_register_test() -> None:
    counter = Counter()
    sequences = ['AaTg', 'aAtG']
    counts_to_test = counter.count_in_sequences(sequences)
    sample_counts = {
        0: {'A': 2, 'T': 0, 'G': 0},
        1: {'A': 2, 'T': 0, 'G': 0},
        2: {'A': 0, 'T': 2, 'G': 0},
        3: {'A': 0, 'T': 0, 'G': 2},
    }
    assert counts_to_test == sample_counts, 'нихера не работает, переделывай'


def hyphen_test() -> None:
    counter = Counter()
    sequences = ['AA-AA-GG', '--A-G-']
    counts_to_test = counter.count_in_sequences(sequences)
    sample_counts = {
        0: {'A': 1, 'G': 0},
        1: {'A': 1, 'G': 0},
        2: {'A': 1, 'G': 0},
        3: {'A': 1, 'G': 0},
        4: {'A': 1, 'G': 1},
        5: {'A': 0, 'G': 0},
        6: {'A': 0, 'G': 1},
        7: {'A': 0, 'G': 1},
    }
    assert counts_to_test == sample_counts, 'нихера не работает, переделывай'
