import re


class Counter:
    """Class for counting unique nucleotides in DNA and RNA sequences."""

    def __init__(self, **kwargs):
        self.vocabulary = kwargs.get('vocabulary')
        self.count_dict = None

    @staticmethod
    def get_vocabulary(sequences: list[str]) -> list[str]:
        """Get unique tokens from list of sequences
        Args:
            sequences: list of sequences
        Returns:
            list of unique tokens in sequences
        """
        uniq_tokens: set[str] = set()
        for seq in sequences:
            seq_tokens = set(seq.upper())
            uniq_tokens.update(seq_tokens)

        if '-' in uniq_tokens:
            uniq_tokens.remove('-')

        return list(uniq_tokens)

    @staticmethod
    def get_sorted_sequences(sequences: list[str]) -> list[str]:
        """Sort sequences according to their length in descending order
        Args:
            sequences: list of sequences
        Returns:
            sorted list of sequences
        """
        return sorted(sequences, key=lambda x: len(x), reverse=True)

    @staticmethod
    def clean_sequences(sequences: list[str]) -> list[str]:
        """Clean sequences from random characters using regex
        Args:
            sequences: list of sequences
        Returns:
            list of cleaned sequences
        """
        clean_sequences: list[str] = [
            ''.join(re.findall(r'[a-zA-Z]+|[-]+', sequence))
            for sequence in sequences
        ]
        return clean_sequences

    def count_one_sequence(self, sequence: str) -> None:
        """Count token occuranse in one sequense and store in class attribute
        Args:
            sequence: sequence to be analyzed
        Returns:
            None
        """
        split_sequence = list(sequence.upper())
        for idx in range(len(split_sequence)):
            if split_sequence[idx] != '-':
                self.count_dict[idx][split_sequence[idx]] += 1

    def get_max_len(self, sequences: list[str]) -> int:
        """Get maximum length of sequences in input list
        Args:
            sequences: list of sequences
        Returns:
            maximum length of sequences
        """

        sorted_sequences = self.get_sorted_sequences(sequences)
        return len(sorted_sequences[0])

    def get_count_dict(
        self,
        sequences: list[str]
    ) -> dict[int, dict[str, int]]:
        """Get dictionary to count token occurance
        Args:
            sequences: list of sequences
        Returns:
            dictionary for encouning tokens
        """
        positions = range(self.get_max_len(sequences))
        count_dict: dict[int, dict[str, int]] = {
            position: dict(
                zip(self.vocabulary, [0] * len(self.vocabulary))
            ) for position in positions
        }
        return count_dict

    def count_in_sequences(
        self,
        sequences: list[str]
    ) -> dict[int, dict[str, int]]:
        """Count token occurance in multiple sequences
        Args:
            sequences: list of sequences
        Returns:
            dictionary containing token counts
        """
        sequences: list[str] = self.clean_sequences(sequences)

        if self.vocabulary is None:
            self.vocabulary = self.get_vocabulary(sequences)

        self.count_dict: dict[int, dict[str, int]] = self.get_count_dict(
            sequences
        )

        for sequence in sequences:
            self.count_one_sequence(sequence)

        return self.count_dict.copy()
