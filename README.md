# Py-Brisk

Py-Brisk is the Python implementation of Brisk, a C++ library built for k-mer indexing.

## CLI

```
# Create an index with k=3 and show statistics
python kmer_cli.py -k 3 index input.fasta --stats

# Create an index and compute minimizers with window size 10
python kmer_cli.py -k 3 index input.fasta --minimizers 10

# Query specific k-mers
python kmer_cli.py -k 3 query ATC GTA CGT

# Use compact index mode
python kmer_cli.py -k 3 --compact index input.fasta
```