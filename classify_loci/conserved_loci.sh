#! /bin/bash

/home/software/Python3/bin/python3 ./conserved_loci.py --species_split Acipenser_ruthenus,Polyodon_spathula \
                                                    -out_folder conserved_loci \
                                                    --species_position  1


exit 0

#--Note--#

#1) input folder name (--species_split A,B,C,D,E)
#2) folder name of output (-out_folder common_loci)
#3) which folder to generate fasta (--species_position 1) means A 
