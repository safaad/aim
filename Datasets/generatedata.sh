#! /usr/bin/bash 

set -x

seqLengths=("100" "1000" "5000" "10000" "20000" "25000" "50000" "100000")

for seqLength in ${seqLengths[@]}; do
    ./bin/generate_dataset --n 100000 -l $seqLength --e 0.01 --o ./const-seq/seq-l$seqLength-e1-100KPairs
done

for seqLength in ${seqLengths[@]}; do
    ./bin/generate_dataset --n 100000 -l $seqLength --e 0.05 --o ./const-seq/seq-l$seqLength-e5-100KPairs
done

echo "finished generating data"
