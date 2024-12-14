import unittest
from kmer_index import KmerIndex, CompactKmerIndex


class TestKmerIndex(unittest.TestCase):
    def setUp(self):
        self.index = KmerIndex(k=3)
        self.test_seq = "ATCGTA"
        self.test_seq_kmers = {"ATC", "TCG", "CGT", "GTA"}

    def test_kmer_generation(self):
        """Test that k-mers are correctly generated from a sequence."""
        generated_kmers = {
            kmer for kmer, _ in self.index._generate_kmers(self.test_seq)
        }
        self.assertEqual(generated_kmers, self.test_seq_kmers)

    def test_sequence_addition(self):
        """Test adding a sequence to the index."""
        self.index.add_sequence(self.test_seq)
        for kmer in self.test_seq_kmers:
            self.assertGreater(len(self.index.query_kmer(kmer)), 0)

    def test_position_tracking(self):
        """Test that k-mer positions are correctly tracked."""
        self.index.add_sequence(self.test_seq)
        # ATC should be at position 0
        self.assertTrue((0, 0) in self.index.query_kmer("ATC"))
        # GTA should be at position 3
        self.assertTrue((0, 3) in self.index.query_kmer("GTA"))

    def test_case_insensitivity(self):
        """Test that k-mer queries are case-insensitive."""
        self.index.add_sequence(self.test_seq)
        self.assertEqual(self.index.query_kmer("ATC"), self.index.query_kmer("atc"))

    def test_multiple_sequences(self):
        """Test handling of multiple sequences."""
        seq1 = "ATCGTA"
        seq2 = "ATCGTA"  # Identical sequence
        self.index.add_sequence(seq1, sequence_id=0)
        self.index.add_sequence(seq2, sequence_id=1)

        atc_positions = self.index.query_kmer("ATC")
        self.assertEqual(len(atc_positions), 2)
        self.assertTrue((0, 0) in atc_positions)
        self.assertTrue((1, 0) in atc_positions)

    def test_n_handling(self):
        """Test that k-mers containing N are properly handled."""
        seq_with_n = "ATNGTA"
        self.index.add_sequence(seq_with_n)
        # Only GTA should be indexed
        self.assertEqual(len(self.index.query_kmer("ATN")), 0)
        self.assertGreater(len(self.index.query_kmer("GTA")), 0)

    def test_minimizers(self):
        """Test minimizer computation."""
        sequence = "ATCGTATCG"
        window_size = 5
        minimizers = self.index.get_minimizers(sequence, window_size)
        self.assertTrue(len(minimizers) > 0)
        # Check that all minimizers are of length k
        self.assertTrue(all(len(m) == self.index.k for m in minimizers))

    def test_statistics(self):
        """Test statistics computation."""
        self.index.add_sequence(self.test_seq)
        stats = self.index.get_statistics()
        self.assertEqual(stats["k"], 3)
        self.assertEqual(stats["unique_kmers"], len(self.test_seq_kmers))
        self.assertEqual(stats["sequences_indexed"], 1)


class TestCompactKmerIndex(unittest.TestCase):
    def setUp(self):
        self.index = CompactKmerIndex(k=3)
        self.test_seq = "ATCGTA"

    def test_kmer_encoding(self):
        """Test binary encoding/decoding of k-mers."""
        original_kmer = "ATC"
        encoded = self.index._encode_kmer(original_kmer)
        decoded = self.index._decode_kmer(encoded)
        self.assertEqual(original_kmer, decoded)

    def test_sequence_addition(self):
        """Test adding sequence to compact index."""
        self.index.add_sequence(self.test_seq)
        self.assertGreater(len(self.index.kmer_array), 0)
        self.assertEqual(len(self.index.kmer_array), len(self.index.position_array))


if __name__ == "__main__":
    unittest.main()
