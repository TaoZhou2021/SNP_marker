#! /bin/bash

species='albus platorycnhsu platorynchus2 suttkusi'

#1. reads_assembly
for i in $species
do 
  ./reads_assembly $i
done

#2. call_consensus
./call_consensus snp 2

#3. snp_calling
./call_snp snp

exit 0
