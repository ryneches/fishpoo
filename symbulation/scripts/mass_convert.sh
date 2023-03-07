#! /usr/bin/bash

for ((i=${1}; i<=${2}; i++))
do
    echo ./symbulation_to_abstract_newick_pipeline.sh ../alife-std-dev-python/data/ $i $3 $4
    ./symbulation_to_abstract_newick_pipeline.sh ../alife-std-dev-python/data/ $i $3 $4
done