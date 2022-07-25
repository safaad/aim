We evaluate AIM using two types of datasets: (1) real datasets for short sequence pairs, (2) synthetic datasets for long sequence pairs

## Real datasets for short sequence pairs (100bp, 150bp, and 250bp long)
The short read-reference pair sets are generated using [minimap2](https://lh3.github.io/minimap2/) by mapping the following datasets to the human reference genome GRCh37:

  1. https://www.ebi.ac.uk/ena/browser/view/ERR240727
  2. https://www.ebi.ac.uk/ena/browser/view/SRR826460
  3. https://www.ebi.ac.uk/ena/browser/view/SRR826471

The human reference genome can be downloaded from:
```bash
ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
```

After extracting the read-reference pairs we format the datasets as shown in the samples above.

## Synthetic datasets for long sequence pairs (500bp, 1000bp, 5Kbp, and 10Kbp long)
We use the data generator found in WFA's repository (https://github.com/smarco/WFA) to simulate long sequence pairs:
```bash
git clone https://github.com/smarco/WFA
cd WFA
make
# Example of generating synthetic dataset with sequence length 100, % edit distance 1% and 5M pairs
./bin/generate_dataset --n 5000000 --l 100 --e 0.01 --o ./synthetic-l100-e1-5MPairs
```