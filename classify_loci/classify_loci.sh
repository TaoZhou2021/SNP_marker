#! /bin/bash

/home/software/Python3/bin/python3 ./classify_loci.py --input_file1 Lepisosteus_oculatus.dna.fas \
                                                      -out_folder classify_loci \
                                                      -rf  Acipenser_ruthenus-60scaffolds.fa \
                                                      -output1 single-copy_loci.fa \
                                                      -output2 double-copy_loci.fa \
                                                      -e 0.0001 \
                                                      -db Acipenser \
                                                      -bf Lep2Acipenser-blast.out                               


exit 0

#--Note--#

#1) input baits fasta file (-inputfile1 Lepisosteus_oculatus.dna.fas)
#2) folder name of output (-out_folder classify_loci)
#3) reference genome file (-rf Acipenser_ruthenus-60scaffolds.fa)
#4) database prefix (-db Acipenser)
#5) blastn parameter e value (-evalue 0.0001)
#6) blastn result (-bf Lep2Acipenser-blast.out)
#7) single-copy loci fasta (-output1 single-copy_loci.fa)
#8) double-copy loci fasta (-output2 double-copy_loci.fa)
