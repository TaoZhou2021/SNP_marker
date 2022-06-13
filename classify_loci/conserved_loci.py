#! /home/software/Python3/bin/python3
# coding: utf-8

import os,sys
import argparse
import pysam
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = argparse.ArgumentParser(description='align sequences in two fasta files')
parser.add_argument('--input_file1', '-input1', dest='inputfile1', help='input file1')
parser.add_argument('--input_file2', '-input2', dest='inputfile2', help='input file2')
parser.add_argument('--output_file', '-output', dest='outputfile', help='output fasta file')
parser.add_argument('--classify_file', '-cf', dest='classifyfile', help='classify file')
parser.add_argument('--output_file1', '-output1', dest='outputfile1', help='output txt file1')
parser.add_argument('--output_file2', '-output2', dest='outputfile2', help='output txt file2')
parser.add_argument('--output_file3', '-output3', dest='outputfile3', help='output txt file3')
parser.add_argument('--output_file4', '-output4', dest='outputfile4', help='output txt file4')
parser.add_argument('--database_prefix', '-db', dest='db_prefix', help='prefix of database')
parser.add_argument('--folder_names', '-f', dest='folder_names', help='names of folder')
parser.add_argument('--species_split','-ss', dest='species_split', help='split of species name')
parser.add_argument('--species_position','-sp', dest='species_position', help='position of species name')
parser.add_argument('--loci_type', '-loci', dest='loci_type', help='type of loci')
parser.add_argument('--E_value', '-evalue', dest='e_value', help='thread of e value')
parser.add_argument('--blast_file', '-bf', dest='blastn_output', help='output of blastn result')
parser.add_argument('--ref_file', '-rf', dest='ref_file', help='reference file')
parser.add_argument('--output_floder', '-out_folder', dest='outputfolder', help='folder of the output') 

args = parser.parse_args()

folder0='/3-loci/'
folder1='/1-conserved_loci-txt/'
folder2='/2-conserved_loci-fa/'

dict1={}
dict_chr={}

def make_dir(file_path):
    if os.path.exists(file_path):
        pass
    else:
        os.makedirs(file_path)

def conserved_loci(species_split,loci_type,outputfolder):
    outfolder=outputfolder+folder1
    make_dir(outfolder)
    split=species_split.split(',')
    num=len(split)
    common_list,common_dict=[],{}
    file_dict={'single-copy':'single-copy_loci.fa','double-copy':'double-copy_loci.fa'}
    add_dict= {'single-copy':1, 'double-copy':0.5}       
    for i in range(0,num):
        inputfile=split[i] + folder0 + file_dict[loci_type]         
        with open(inputfile,'r') as f1:
            for title,seq in SimpleFastaParser(f1):
                loci=title.split()[0].split('|')[0]
                if loci not in common_dict.keys():
                    common_dict[loci]=add_dict[loci_type]
                else:
                    common_dict[loci]+=add_dict[loci_type]
    outputfile=outfolder+'common-'+loci_type+'_loci.txt' 
    with open(outputfile,'w') as f2: 
        for key,value in common_dict.items():
            if value==num:
                f2.writelines("%s\n" % (key))
                common_list.append(key)
    return common_list

def pick_loci(species_split,species_position,outputfolder):
    outfolder= outputfolder+folder2
    make_dir(outfolder)
    pos=int(species_position) -1
    split=species_split.split(',')
    species=split[pos]
    loci_name=['single-copy','double-copy']
    for i in loci_name: 
        file_list=conserved_loci(species_split,i,outputfolder)
        input_file=species + folder0 + i + '_loci.fa'
        output_file=outfolder + species + '-conserved-' + i + '_loci.fa'
        with open(input_file,'r') as f1,open(output_file,'w') as f2:
            for title,seq in SimpleFastaParser(f1):
                loci = title.split()[0].split('|')[0]
                if loci in file_list:
                    f2.writelines(">%s\n%s\n" % (title,seq))  

if __name__=="__main__": 
 if "--input_file1" or "-input1" in str(sys.argv):
  pick_loci(args.species_split,args.species_position,args.outputfolder)
