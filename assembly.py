#! /home/software/Python3/bin/python3
# coding: utf-8

import sys,os
import argparse
import itertools
import shutil
import subprocess
from decimal import Decimal
from Bio.Seq import Seq
from Bio import SeqIO,Align
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


parser = argparse.ArgumentParser(description='align sequences in two fasta files')
parser.add_argument('--input_fq1', '-inputfq1', dest='inputfq1', help='input fastq file1')
parser.add_argument('--input_fq2', '-inputfq2', dest='inputfq2', help='input fastq file2')
parser.add_argument('--database_prefix', '-db', dest='db_prefix', help='prefix of database')
parser.add_argument('--E_value', '-evalue', dest='e_value', help='thread of e value')
parser.add_argument('--blast_file', '-bf', dest='blastn_output', help='output of blastn result')
parser.add_argument('--ref_file', '-rf', dest='ref_file', help='reference file')
parser.add_argument('--ref_file1', '-rf1', dest='ref_file1', help='reference file1')
parser.add_argument('--ref_file2', '-rf2', dest='ref_file2', help='reference file2')
parser.add_argument('--output_fa1', '-outputfa1', dest='outputfa1', help='output fasta file1')
parser.add_argument('--output_fa2', '-outputfa2', dest='outputfa2', help='output fasta file2')
parser.add_argument('--output_floder', '-out_folder', dest='outputfolder', help='folder of the output')
parser.add_argument('--g1_path', '-g1_folder', dest='g1folder', help='path of g1 folder')
parser.add_argument('--g2_path', '-g2_folder', dest='g2folder', help='path of g2 folder')
args = parser.parse_args()

folder1='/1-oneline/'
folder2='/2-blast/'
folder3='/3-pickfq/'
folder4='/4-assemble/'

def make_dir(file_path):
    if os.path.exists(file_path):
        pass
    else:
        os.makedirs(file_path)

def fq2fa(inputfile1, inputfile2, outputfolder): 
    fa_folder = outputfolder + folder1
    make_dir(fa_folder)
    outputfile1 = fa_folder + os.path.splitext(os.path.basename(inputfile1))[0]+'.fa'
    outputfile2 = fa_folder + os.path.splitext(os.path.basename(inputfile2))[0]+'.fa'
    SeqIO.convert(inputfile1, "fastq", outputfile1, "fasta")
    SeqIO.convert(inputfile2, "fastq", outputfile2, "fasta")
    print("convert fastq to fasta : OK")
    return outputfile1, outputfile2

def one_line(inputfq1,inputfq2,outputfolder):
    inputfa1, inputfa2 = fq2fa(inputfq1, inputfq2, outputfolder)
    seq_dict1={}
    seq_dict2={}
    outputfile = outputfolder + folder1 + os.path.splitext(os.path.basename(inputfq1))[0][:-3]+'-oneline-merge.fa'    
    with open(inputfa1,'r') as f1:
        with open(inputfa2,'r') as f2:     
            with open(outputfile,'w') as f:
                lines1=f1.readlines()
                lines2=f2.readlines()
                for line1 in lines1:
                    if line1.startswith('>'):
                        name1=line1.split()[0]
                        seq_dict1[name1]=''
                    else:
                        seq_dict1[name1]+=line1.replace('\n','')
                for line2 in lines2:
                    if line2.startswith('>'):
                        name2=line2.split()[0]
                        seq_dict2[name2]=''
                    else:
                        seq_dict2[name2]+=line2.replace('\n','')      
                for key,value in seq_dict1.items():
                    seq2=seq_dict2[key]
                    newseq=value+'NNNNN'+ Seq(seq2).reverse_complement()
                    f.writelines("%s\n%s\n" % (key,newseq))
    print("merge fasta to one line : OK")    

def rmdup(inputfq1,outputfolder):
    inputfa = outputfolder + folder1 + os.path.splitext(os.path.basename(inputfq1))[0][:-3]+'-oneline-merge.fa' 
    outputfa = outputfolder + folder1 + os.path.splitext(os.path.basename(inputfq1))[0][:-3]+'-oneline-merge-rmdup.fa'
    rmdup_cline = "seqkit rmdup " + inputfa + " -s -o " + outputfa
    subp = subprocess.Popen(str(rmdup_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    subp.wait()
    if subp.poll() == 0:
        print("rmdup OK")
    else:
        print("rmdup failed")                

def make_db(ref_file,db_prefix,outputfolder):
    blast_folder=outputfolder + folder2
    make_dir(blast_folder)
    makedb_cline="makeblastdb -in " + ref_file + " -dbtype nucl -parse_seqids -out " + blast_folder + db_prefix
    subp = subprocess.Popen(str(makedb_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    subp.wait()
    if subp.poll() == 0:
        print("makedb OK")
    else:
        print("makedb failed")

def blastn2classify_reads(inputfq1,inputfq2,outputfolder,ref_file,db_prefix,e_value,blastn_output):
    inputfa = outputfolder + folder1 + os.path.splitext(os.path.basename(inputfq1))[0][:-3]+'-oneline-merge-rmdup.fa'
    inputdb = outputfolder + folder2 + db_prefix
    blastn_output_path = outputfolder + folder2
    cmd_path = '/home/room/users/zhoutao/ncbi-blast-2.11.0+/bin/blastn'    
    blastn_cline = cmd_path + " -task blastn -query " + inputfa + " -out " + blastn_output_path + blastn_output + " -db " + inputdb + " -outfmt 6 -num_threads 30 -evalue " + e_value
    subp = subprocess.Popen(str(blastn_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    subp.wait()
    if subp.poll() == 0:
        print("blastn OK")
    else:
        print("blastn failed")       
    classify_output=blastn_output_path + os.path.splitext(os.path.basename(inputfa))[0] + '_reads2loci.txt' 
    with open(blastn_output_path+blastn_output,'r') as f1:
        with open(classify_output,'w') as f:
            lines=f1.readlines()
            read_dict={}
            for line in lines:
                split=line.split('\t')
                read=split[0]
                loci=split[1]
                score=int(split[3]) - 2*int(split[4]) - 2.5*int(split[5])
                if read not in read_dict.keys():
                    read_dict[read]=[[loci,score]]
                    loci1=loci
                    score1=score
                else:
                    if score > score1 :
                        read_list = read_dict[read]
                        read_list.insert(0,[loci,score])
                        score1=score
                    elif score == score1:
                        read_list = read_dict[read]
                        read_list.insert(1,[loci,score])      
            s=sorted(read_dict.items(),key=lambda x:x[0])    
            for key,value in s:
                if len(value)== 1:
                    loci1=value[0][0]
                    f.writelines("%s\t%s\n" % (key,loci1))
                else:
                    loci1=value[0][0]
                    score1=int(value[0][1])
                    loci2=value[1][0]
                    score2=int(value[1][1])
                    if score1 == score2:
                        f.writelines("%s\t%s\n" % (key,loci1))
                        f.writelines("%s\t%s\n" % (key,loci2)) 
    print("classify reads : OK")      
    return classify_output

def file_name_listdir_local(file_dir):
    files_local = []
    for files in os.listdir(file_dir):
        files_local.append(files)
    return files_local

def txt2list(inputfile):
    with open(inputfile,'r') as f:
        loci_list=[]
        lines=f.readlines()
        for line in lines:
            name=line.split()[0]
            if name !='' and name !='\n':
                if name not in loci_list:
                    loci_list.append(name)
    return loci_list

def make_txt(inputfile,outputfolder):
    folder_path = outputfolder + folder3 + 'loci-txt_fq/'
    make_dir(folder_path)
    with open(inputfile,'r') as f1:
        lines=f1.readlines()
        dict1 = {}
        list1 = []
        for line in lines:
            split=line.split()
            seq1=split[0]
            loci1=split[1]
            if loci1 not in list1:
                list1.append(loci1)
                dict1[loci1]=[seq1]
            else:
                dict1[loci1].append(seq1)
        for key,value in dict1.items():
            txt_name= key + '.txt'
            f3=open(folder_path+txt_name,'w')
            for i in value:
                f3.writelines("%s\n" % (i))
            f3.close()
    return folder_path 

def make_fq(inputfq1,inputfq2,inputfile,outputfolder):
    folder_path = make_txt(inputfile,outputfolder)
    with open(inputfq1,'r') as fq1:
        with open(inputfq2,'r') as fq2:
            seq1_dict=dict((title.split()[0],[seq,qual]) for title,seq,qual in FastqGeneralIterator(fq1))      
            seq2_dict=dict((title.split()[0],[seq,qual]) for title,seq,qual in FastqGeneralIterator(fq2)) 
            txt_list=file_name_listdir_local(outputfolder + folder3 + 'loci-txt_fq')
            for i in txt_list:
                if i.endswith('.txt'):
                    fq_list=txt2list(folder_path+i)
                    prefix=i[:-4]
                    fq1_name=folder_path+prefix+'.txt_R1.fq' 
                    fq2_name=folder_path+prefix+'.txt_R2.fq'
                    with open(fq1_name,'w') as F1: 
                        with open(fq2_name,'w') as F2:
                            for j in fq_list: 
                                seq_R1=seq1_dict[j][0]
                                qual_R1=seq1_dict[j][1]
                                F1.writelines("@%s\n%s\n+\n%s\n" % (j,seq_R1,qual_R1))
                                seq_R2=seq2_dict[j][0]
                                qual_R2=seq2_dict[j][1]
                                F2.writelines("@%s\n%s\n+\n%s\n" % (j,seq_R2,qual_R2))

def ref_loci(inputfa1):
    with open(inputfa1,'r') as f1:
         ref_list=[]
         lines=f1.readlines()
         for line in lines:
             if line.startswith(">"):
                 name=line.split()[0][1:]
                 if name not in ref_list:
                     ref_list.append(name)
         return ref_list

def cp_file(inputfa1,fq_folder,outputfolder):
    file_list=ref_loci(inputfa1)
    folder_path = outputfolder + folder3 + 'loci-txt_fq/'
    output_path = outputfolder + folder3 + fq_folder + '/'
    make_dir(output_path)
    for i in file_list:
        fq1_name=folder_path+i+'.txt_R1.fq'
        fq2_name=folder_path+i+'.txt_R2.fq'
        if os.path.exists(fq1_name):
            shutil.copy(fq1_name,output_path)
            shutil.copy(fq2_name,output_path)
        else:
            print("%s not found\n%s not found" % (fq1_name,fq2_name))
            pass

def fq1_name_listdir_local(file_dir):
    files_local = []
    for filename in os.listdir(file_dir):
        if "|" in str(filename):
            newname = filename.replace("|","-")
            os.rename(file_dir+filename, file_dir+newname)
    for files in os.listdir(file_dir):
        name = os.path.splitext(files)[0][-3:]
        if name == '_R1':
            files_local.append(files)
    return files_local

def files_name_listdir_local(file_dir,suffix):
    files_local = []
    for files in os.listdir(file_dir):
        name = os.path.splitext(files)[1][1:]
        size = os.path.getsize(file_dir + files)
        if name == suffix and size != 0:
            files_local.append(files)
    return files_local


def assemble(g_folder,outputfolder):
    input_dir = outputfolder + folder3 + g_folder +'/'
    output_dir = outputfolder + folder4 + g_folder + '/AByss/'
    make_dir(output_dir)
    fq1_list = fq1_name_listdir_local(input_dir)
    cmd_path = '/home/software/abyss-1.3.7/ABYSS/ABYSS'
    for fq1 in fq1_list:
        inputfq1 = input_dir + os.path.basename(fq1)
        inputfq2 = input_dir + os.path.splitext(os.path.basename(fq1))[0][:-3]+'_R2.fq'
        outputfa = output_dir + os.path.splitext(os.path.basename(fq1))[0][:-3]+'.contig'  
        assemble_cline = cmd_path + " -k 31 " + inputfq1 + " "+ inputfq2 + " -o " + outputfa
        subp = subprocess.Popen(str(assemble_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
        subp.wait()    
        if subp.poll ==0:
            print("AByss assemble OK")
        else:
            print("AByss assemble failed")


def score(seqA,seqB):
    aligner = Align.PairwiseAligner(match_score=1.0)
    seq1=str(seqA)
    seq2=str(seqB)
    score = aligner.score(seq1, seq2)
    return score

def best_contig(ref_file,g_folder,species,outputfolder):
    input_dir = outputfolder + folder4 + g_folder + '/AByss/'
    output_dir = outputfolder + folder4 + g_folder + '/score/'
    make_dir(output_dir)
    suffix = 'contig'
    inputfiles = files_name_listdir_local(input_dir,suffix)
    for inputfile1 in inputfiles:
        outputfile1 = output_dir + os.path.splitext(inputfile1)[0] + '.best_contig'
        with open(input_dir + inputfile1,'r') as f1: 
            with open(ref_file,'r') as f2:
                with open(outputfile1,'w') as f3:
                    fa_dict=dict((title,seq) for title,seq in SimpleFastaParser(f2))
                    ref_name=os.path.splitext(os.path.splitext(inputfile1)[0])[0].replace("-","|")
                    ref_fa=fa_dict[ref_name]
                    len_fa=len(ref_fa)
                    lines=f1.readlines()
                    list1 =[]
                    dict1= {}
                    for line in lines:
                        line=line.strip()
                        if line.startswith('>'):
                            name=line
                            dict1[name]=''
                        else:
                            dict1[name]+=line
                    for key,value in dict1.items():
                        score1=score(ref_fa,value)  
                        if list1 ==[]:
                            list1.append([key,score1])
                            score2=score1
                        else:
                            if int(score1) > int(score2) and int(score1) > len_fa * 0.5 :
                                list1.insert(0,[key,score1])
                                score2=score1
                            elif int(score1) == int(score2) and  int(score1) > len_fa * 0.5 :
                                list1.insert(1,[key,score1])
                                score2=score1
                    if len(list1) == 1:
                        score0=list1[0][1]
                        name=list1[0][0]
                        seq=dict1[name]
                        f3.writelines(">contig1 (%s/%s) %s\n%s" % (score0,str(len_fa),species,seq))
                    elif len(list1) >= 2:
                        score2=list1[1][1]
                        if int(score2) > len_fa * 0.5 :
                            score0=list1[0][1]
                            name1=list1[0][0]
                            seq1=dict1[name1]
                            name2=list1[1][0]
                            seq2=dict1[name2] 
                            f3.writelines(">contig1 (%s/%s) %s\n%s\n>contig2 (%s/%s) %s\n%s" % (score0,str(len_fa),species,seq1,score2,str(len_fa),species,seq2))

def make_db2(ref_file,db_prefix,g_folder,outputfolder):
    blast_folder=outputfolder + folder4 + g_folder + '/blast_db/'
    make_dir(blast_folder)
    makedb_cline="makeblastdb -in " + ref_file + " -dbtype nucl -parse_seqids -out " + blast_folder + db_prefix
    subp = subprocess.Popen(str(makedb_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    subp.wait()
    if subp.poll() == 0:
        print("makedb OK")
    else:
        print("makedb failed")

def reblast(db_prefix,e_value,g_folder,outputfolder):
    input_dir = outputfolder + folder4 + g_folder + '/score/'
    output_dir = outputfolder + folder4 + g_folder + '/reblast/'
    make_dir(output_dir)
    files_list = files_name_listdir_local(input_dir,'best_contig')
    inputdb = outputfolder + folder4 + g_folder + '/blast_db/' + db_prefix
    cmd_path = '/home/software/conda/miniconda/bin/blastn'
    for files in files_list:
        inputfa= input_dir + files
        blastn_output = output_dir + os.path.splitext(files)[0] + '.reblast'                                                  
        blastn_cline = cmd_path + " -task blastn -query " + inputfa + " -out " + blastn_output + " -db " + inputdb + " -outfmt 6 -num_threads 30 -evalue " + e_value
        subp = subprocess.Popen(str(blastn_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
        subp.wait()
        if subp.poll() == 0:
            print("reblast OK")
        else:
            print("reblast failed")

def final_contig(species,g_folder,outputfolder):
    input_dir = outputfolder + folder4 + g_folder
    contig_dir = input_dir + '/score/'
    reblast_dir = input_dir + '/reblast/'
    inputfiles1 = files_name_listdir_local(reblast_dir,'reblast')
    output_dir = input_dir + '/final_contig/'
    make_dir(output_dir)
    for inputfile1 in inputfiles1:
        inputfile2 = contig_dir + os.path.splitext(inputfile1)[0]+'.best_contig'
        outputfile = output_dir + os.path.splitext(inputfile1)[0]+'.final_contig'
        with open(inputfile2,'r') as f1:
            with open(reblast_dir + inputfile1,'r') as f2:
                with open(outputfile,'w') as f3:
                    list1=[]
                    list2=[]
                    lines=f2.readlines()
                    for line in lines:
                       if line !='':
                           split=line.split()
                           loci=split[1]
                           if str(loci) + '.txt.reblast' == os.path.basename(inputfile1).replace("-","|"):
                               contig=split[0]
                               match=int(split[3])
                               mismatch=int(split[4])
                               gapopen=int(split[5])
                               start=int(split[6])
                               end=int(split[7])
                               score = match - mismatch - gapopen     
                               list1.append([contig,score,start,end])
                    contig_dict=dict((title.split()[0],seq) for title,seq in SimpleFastaParser(f1))
                    if list1 !=[]:
                        if len(list1) == 1:
                            contig=list1[0][0]
                            start=list1[0][2]
                            end=list1[0][3]
                            seq=contig_dict[contig][start:end]
                            f3.writelines(">%s\t[%s:%s] %s\n%s\n" % (contig,start,end,species,seq))
                        else:
                            score1=list1[0][1]
                            score2=list1[1][1]
                            if score1 > score2 :
                                contig=list1[0][0]
                                start=list1[0][2]
                                end=list1[0][3]         
                                f3.writelines(">%s\t[%s:%s] %s\n%s\n" % (contig,start,end,species,seq))
                            else:
                                contig=list1[1][0]
                                start=list1[1][2]
                                end=list1[1][3]
                                f3.writelines(">%s\t[%s:%s] %s\n%s\n" % (contig,start,end,species,seq))




if __name__=="__main__": 
    if "--input_fq1" or "inputfq1" in str(sys.argv):
        one_line(args.inputfq1,args.inputfq2,args.outputfolder)
        rmdup(args.inputfq1,args.outputfolder)
        make_db(args.ref_file,args.db_prefix,args.outputfolder)
        classify_txt = blastn2classify_reads(args.inputfq1,args.inputfq2,args.outputfolder,args.ref_file,args.db_prefix,args.e_value,args.blastn_output)
        make_fq(args.inputfq1,args.inputfq2,classify_txt,args.outputfolder)
        cp_file(args.ref_file1,args.g1folder,args.outputfolder)
        cp_file(args.ref_file2,args.g2folder,args.outputfolder)
        assemble(args.g1folder,args.outputfolder)
        assemble(args.g2folder,args.outputfolder)
        best_contig(args.ref_file1,args.g1folder,args.speciesname,args.outputfolder)
        best_contig(args.ref_file2,args.g2folder,args.speciesname,args.outputfolder)
        make_db2(args.ref_file1,args.db_prefix1,args.g1folder,args.outputfolder)
        make_db2(args.ref_file2,args.db_prefix2,args.g2folder,args.outputfolder)
        reblast(args.db_prefix1,args.e_value,args.g1folder,args.outputfolder)
        reblast(args.db_prefix2,args.e_value,args.g2folder,args.outputfolder)
        final_contig(args.speciesname,args.g1folder,args.outputfolder)
        final_contig(args.speciesname,args.g2folder,args.outputfolder)
