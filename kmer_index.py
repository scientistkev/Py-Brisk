from collections import defaultdict
import itertools
from typing import Dict, Set, Iterator, List
import numpy as np


class KmerIndex:
    def __init__(self, k: int):
        """Initialize a k-mer index with specified k-mer length.

        Args:
            k (int): Length of k-mers to index
        """
        self.k = k
        self.kmer_dict = defaultdict(set)  # k-mer to position mapping
        self.sequence_count = 0

    def _generate_kmers(self, sequence: str) -> Iterator[tuple[str, int]]:
        """Generate k-mers from a sequence along with their positions.

        Args:
            sequence (str): Input DNA sequence

        Yields:
            tuple[str, int]: (k-mer, position) pairs
        """
        sequence = sequence.upper()
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i : i + self.k]
            if "N" not in kmer:  # Skip k-mers containing N
                yield kmer, i

    def add_sequence(self, sequence: str, sequence_id: int = None) -> None:
        """Add a sequence to the index.

        Args:
            sequence (str): DNA sequence to index
            sequence_id (int, optional): ID for the sequence. If None, auto-increment
        """
        if sequence_id is None:
            sequence_id = self.sequence_count
            self.sequence_count += 1

        for kmer, pos in self._generate_kmers(sequence):
            self.kmer_dict[kmer].add((sequence_id, pos))

    def query_kmer(self, kmer: str) -> Set[tuple[int, int]]:
        """Query the index for a specific k-mer.

        Args:
            kmer (str): K-mer to query

        Returns:
            Set[tuple[int, int]]: Set of (sequence_id, position) pairs
        """
        return self.kmer_dict[kmer.upper()]

    def get_minimizers(self, sequence: str, window_size: int) -> List[str]:
        """Get minimizers for a sequence (useful for reducing index size).

        Args:
            sequence (str): Input sequence
            window_size (int): Size of the window for minimizer selection

        Returns:
            List[str]: List of minimizer k-mers
        """
        minimizers = []
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i : i + window_size]
            kmers = list(self._generate_kmers(window))
            if kmers:  # If any valid k-mers in window
                min_kmer = min(kmers, key=lambda x: x[0])[0]
                minimizers.append(min_kmer)
        return minimizers

    def get_statistics(self) -> Dict:
        """Get statistics about the index.

        Returns:
            Dict: Dictionary containing index statistics
        """
        num_kmers = len(self.kmer_dict)
        total_positions = sum(len(positions) for positions in self.kmer_dict.values())
        return {
            "k": self.k,
            "unique_kmers": num_kmers,
            "total_positions": total_positions,
            "sequences_indexed": self.sequence_count,
            "average_positions_per_kmer": (
                total_positions / num_kmers if num_kmers > 0 else 0
            ),
        }


class CompactKmerIndex:
    """Memory-efficient version using bit-packed k-mers."""

    def __init__(self, k: int):
        self.k = k
        self._encoding = {"A": 0, "C": 1, "G": 2, "T": 3}
        self._decoding = {v: k for k, v in self._encoding.items()}
        self.kmer_array = np.zeros((0,), dtype=np.uint64)
        self.position_array = np.zeros((0,), dtype=np.uint32)

    def _encode_kmer(self, kmer: str) -> int:
        """Convert k-mer string to binary representation."""
        encoded = 0
        for base in kmer.upper():
            encoded = (encoded << 2) | self._encoding[base]
        return encoded

    def _decode_kmer(self, encoded: int) -> str:
        """Convert binary representation back to k-mer string."""
        kmer = []
        mask = 3  # Binary mask: 11
        for _ in range(self.k):
            kmer.append(self._decoding[encoded & mask])
            encoded >>= 2
        return "".join(reversed(kmer))

    def add_sequence(self, sequence: str) -> None:
        """Add a sequence to the compact index."""
        new_kmers = []
        new_positions = []

        for kmer, pos in self._generate_kmers(sequence):
            encoded_kmer = self._encode_kmer(kmer)
            new_kmers.append(encoded_kmer)
            new_positions.append(pos)

        if new_kmers:
            self.kmer_array = np.append(self.kmer_array, new_kmers)
            self.position_array = np.append(self.position_array, new_positions)
