import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from kmer_index import KmerIndex, CompactKmerIndex


def parse_args():
    parser = argparse.ArgumentParser(description="K-mer indexing tool")
    parser.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        required=True,
        help="Length of k-mers to generate",
    )
    parser.add_argument(
        "--compact", action="store_true", help="Use memory-efficient compact index"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Index command
    index_parser = subparsers.add_parser("index", help="Create index from FASTA file")
    index_parser.add_argument("input_file", type=Path, help="Input FASTA file")
    index_parser.add_argument(
        "--minimizers", type=int, help="Window size for minimizer computation"
    )
    index_parser.add_argument(
        "--stats", action="store_true", help="Show indexing statistics"
    )

    # Query command
    query_parser = subparsers.add_parser("query", help="Query k-mers")
    query_parser.add_argument("kmers", nargs="+", help="K-mers to query")

    return parser.parse_args()


def process_fasta(index, fasta_path):
    """Process FASTA file and add sequences to index."""
    try:
        for seq_record in SeqIO.parse(str(fasta_path), "fasta"):
            index.add_sequence(str(seq_record.seq), sequence_id=seq_record.id)
    except Exception as e:
        print(f"Error processing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    args = parse_args()

    # Initialize appropriate index type
    IndexClass = CompactKmerIndex if args.compact else KmerIndex
    index = IndexClass(k=args.kmer_size)

    if args.command == "index":
        # Process input file
        process_fasta(index, args.input_file)

        # Calculate minimizers if requested
        if args.minimizers:
            with open(args.input_file) as f:
                for seq_record in SeqIO.parse(f, "fasta"):
                    minimizers = index.get_minimizers(
                        str(seq_record.seq), args.minimizers
                    )
                    print(f"Minimizers for {seq_record.id}:")
                    for min_kmer in minimizers:
                        print(f"  {min_kmer}")

        # Show statistics if requested
        if args.stats and isinstance(index, KmerIndex):
            stats = index.get_statistics()
            print("\nIndexing Statistics:")
            for key, value in stats.items():
                print(f"{key}: {value}")

    elif args.command == "query":
        # Validate k-mer lengths
        invalid_kmers = [k for k in args.kmers if len(k) != args.kmer_size]
        if invalid_kmers:
            print(
                f"Error: These k-mers are not of length {args.kmer_size}:",
                file=sys.stderr,
            )
            print(", ".join(invalid_kmers), file=sys.stderr)
            sys.exit(1)

        # Perform queries
        for kmer in args.kmers:
            print(f"\nPositions for k-mer {kmer}:")
            if isinstance(index, KmerIndex):
                positions = index.query_kmer(kmer)
                for seq_id, pos in positions:
                    print(f"  Sequence {seq_id}, position {pos}")
            else:
                print("  Query not supported for compact index")


if __name__ == "__main__":
    main()
