#! /usr/bin/bash

directory=$1
seed=$2
update=$3
prefix=_${4}_${update}_${seed}_

python abstract_phylogenies.py ${directory}/${seed}/SymSnapshot_Phylogeny__data_UPDATE${update}_SEED${seed}.data ${directory}/${seed}/HostSnapshot_Phylogeny__data_UPDATE${update}_SEED${seed}.data ${directory}/${seed}/InteractionSnapshot_Phylogeny__data_UPDATE${update}_SEED${seed}.data 1000 $prefix
python make_links.py abstract_interactions.csv ${prefix}links.csv

alifedata-phyloinformatics-convert fromalifedata --input-file ${prefix}compressed_sym.csv --output-file ${prefix}symtree.nwk --output-schema newick --suppress-unifurcations
alifedata-phyloinformatics-convert fromalifedata --input-file ${prefix}compressed_host.csv --output-file ${prefix}hosttree.nwk --output-schema newick --suppress-unifurcations