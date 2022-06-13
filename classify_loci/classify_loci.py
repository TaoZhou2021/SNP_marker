#! /home/software/Python3/bin/python3
# coding: utf-8

import os,sys
import argparse
import pysam
import subprocess

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
parser.add_argument('--E_value', '-evalue', dest='e_value', help='thread of e value')
parser.add_argument('--blast_file', '-bf', dest='blastn_output', help='output of blastn result')
parser.add_argument('--ref_file', '-rf', dest='ref_file', help='reference file')
parser.add_argument('--output_floder', '-out_folder', dest='outputfolder', help='folder of the output') 

args = parser.parse_args()

folder1='/1-blastn_db/'
folder2='/2-blastn_output/'
folder3='/3-loci/'

dict1={}
dict_chr={}

def make_dir(file_path):
    if os.path.exists(file_path):
        pass
    else:
        os.makedirs(file_path)

def make_db(ref_file,db_prefix,outputfolder):
    blast_folder=outputfolder + folder1
    make_dir(blast_folder)
    makedb_cline="makeblastdb -in " + ref_file + " -dbtype nucl -parse_seqids -out " + blast_folder + db_prefix
    subp = subprocess.Popen(str(makedb_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    subp.wait()
    if subp.poll() == 0:
        print("makedb OK")
    else:
        print("makedb failed")


def blastn2classify_reads(inputfa,outputfolder,ref_file,db_prefix,e_value,blastn_output):
    inputdb = outputfolder + folder1 + db_prefix
    blastn_output_path = outputfolder + folder2
    make_dir(blastn_output_path) 
    cmd_path = '/home/room/users/zhoutao/ncbi-blast-2.11.0+/bin/blastn'    
    blastn_cline = cmd_path + " -task blastn -query " + inputfa + " -out " + blastn_output_path + blastn_output + " -db " + inputdb + " -outfmt 6 -num_threads 30 -evalue " + e_value
    subp = subprocess.Popen(str(blastn_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    subp.wait()
    if subp.poll() == 0:
        print("blastn OK")
    else:
        print("blastn failed") 


def get_fasta(blastn_output,ref_file,outputfile1, outputfile2, outputfolder):
    loci2_list=[]
    dict1,single_dict,double_dict={},{},{}
    inputfile=outputfolder + folder2+blastn_output
    outfolder=outputfolder + folder3
    make_dir(outfolder)
    output1 = outputfolder + folder2 + os.path.splitext(blastn_output)[0]+'_single-loci.txt'
    output2 = outputfolder + folder2 + os.path.splitext(blastn_output)[0]+'_double-loci.txt'
    output3 = outfolder + outputfile1
    output4 = outfolder + outputfile2
    with open(inputfile,'r') as f1,open(output1,'w') as f2,open(output2,'w') as f3:
        lines=f1.readlines()
        for line in lines:
            split=line.split()
            loci1=split[0]
            name='-'.join(split[:2])
            if loci1 not in dict1.keys():
                dict1[loci1]=1
                loci2_list.append(name)
            else:
                if name not in loci2_list:
                    dict1[loci1]+=1
                else:
                    dict1[loci1]+=2
        for key,value in dict1.items():
            if value==1:
                single_dict[key]=value
            elif value==2:
                double_dict[key]=value
        for line in lines:
            loci1=line.split()[0]
            if loci1 in single_dict.keys():
                f2.writelines("%s\n" % (line.rstrip('\n')))
            elif loci1 in double_dict.keys():
                f3.writelines("%s\n" % (line.rstrip('\n')))
    with open(output1,'r') as f4, open(output2,'r') as f5, open(output3,'w') as f6, open(output4,'w') as f7:
        fasta=pysam.FastaFile(ref_file)
        lines1=f4.readlines()
        for line in lines1:
            split=line.split()
            chr1=str(split[1])
            start=int(split[8])
            end=int(split[9])
            loci=split[0]
            if start < end:
                title=">"+loci+'|'+chr1+'['+str(start)+':'+str(end)+']'+'forward'                 
                seq=fasta.fetch(chr1,start,end)
            else:
                title=">"+loci+'|'+chr1+'['+str(end)+':'+str(start)+']'+'reverse'
                seq=fasta.fetch(chr1,end,start)    
            f6.write("%s\n%s\n" %(title,seq))
        lines2=f5.readlines()
        for line in lines2:
            split=line.split()
            chr1=str(split[1])
            start=int(split[8])
            end=int(split[9])
            loci=split[0]
            if start < end:
                title=">"+loci+'|'+chr1+'['+str(start)+':'+str(end)+']'+'forward'
                seq=fasta.fetch(chr1,start,end)
            else:
                title=">"+loci+'|'+chr1+'['+str(end)+':'+str(start)+']'+'reverse'
                seq=fasta.fetch(chr1,end,start)
            f7.write("%s\n%s\n" %(title,seq))


if __name__=="__main__": 
 if "--input_file1" or "-input1" in str(sys.argv):
  make_db(args.ref_file,args.db_prefix,args.outputfolder)
  blastn2classify_reads(args.inputfile1,args.outputfolder,args.ref_file,args.db_prefix,args.e_value,args.blastn_output)
  get_fasta(args.blastn_output,args.ref_file,args.outputfile1,args.outputfile2,args.outputfolder)
