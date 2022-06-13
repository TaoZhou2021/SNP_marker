
import os,sys
import argparse
import subprocess
import shutil
from decimal import Decimal
from Bio import Align
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import chain
from collections import Counter

parser = argparse.ArgumentParser(description='align sequences in two fasta files')
parser.add_argument('--input_file', '-input', dest='inputfile', help='input txt file')
parser.add_argument('--output_file', '-output', dest='outputfile', help='output txt file')
parser.add_argument('--input_folder','-in_folder', dest='inputfolder', help='folder of inputfile')
parser.add_argument('--output_folder','-out_folder', dest='outputfolder', help='folder of outputfile')
parser.add_argument('--phy_file1','-pf1', dest='phy_file1', help='phyfile1 of fasta')
parser.add_argument('--phy_file2','-pf2', dest='phy_file2', help='phyfile2 of fasta')
parser.add_argument('--nex_file1','-nf1', dest='nex_file1', help='nexfile1 of fasta')
parser.add_argument('--nex_file2','-nf2', dest='nex_file2', help='nexfile2 of fasta')
parser.add_argument('--species_split','-ss', dest='species_split', help='split of species name')
parser.add_argument('--fa_split','-fs', dest='fa_split', help='split of fasta name')
parser.add_argument('--ref_number','-n', dest='ref_number', help='number of reference')
parser.add_argument('--mafft_number','-m', dest='mafft_number', help='number of species to mafft align')
parser.add_argument('--g1_path', '-g1_folder', dest='g1folder', help='path of g1 folder')
parser.add_argument('--g2_path', '-g2_folder', dest='g2folder', help='path of g2 folder')
parser.add_argument('--g1_ref', '-g1_ref', dest='g1ref', help='reference of g1 folder')
parser.add_argument('--g2_ref', '-g2_ref', dest='g2ref', help='reference of g2 folder')
parser.add_argument('--g1_ref2', '-g1_ref2', dest='g1ref2', help='reference2 of g1 folder')
parser.add_argument('--g2_ref2', '-g2_ref2', dest='g2ref2', help='reference2 of g2 folder')
parser.add_argument('--out_group', '-out_group', dest='outgroup', help='fasta of outgroup')

args = parser.parse_args()

folder3='/3-pickfq/'
folder4='/4-assemble/'
folder5='/1-concatenate/'
folder6='/2-mafft_align/'
folder7='/3-fas/'
folder8='/4-call_snp/'

def make_dir(file_path):
    if os.path.exists(file_path):
        pass
    else:
        os.makedirs(file_path)

def files_name_listdir_local(file_dir,suffix):
    files_local = []
    for files in os.listdir(file_dir):
        name = os.path.splitext(files)[1][1:]
        size = os.path.getsize(file_dir + files)
        if name == suffix and size != 0:
            files_local.append(files)
    return files_local

def fa2dict(filepath,inputfile,suffix,g_folder,inputfolder):
    fa_dict={}
    split=filepath.split(':')
    file_number=len(split)
    for i in range(0,file_number):
        file_path=split[i] + inputfolder +  folder4 + g_folder + '/' + suffix + '/'
        fact=os.path.exists(file_path+inputfile)
        if fact is True:
            with open(file_path+inputfile,'r') as f1:
                lines=f1.readlines()
                for line in lines:
                    if line.startswith('>'):
                        name=line.rstrip('\n')
                        fa_dict[name]=''
                    else:
                        fa_dict[name]+=line.replace('\n','')
    return fa_dict

def fa2dict2(inputfile):
    dict1={}
    with open(inputfile,'r') as f0:
        for title,seq in SimpleFastaParser(f0):
            fa_loci=title.split('|')[0]
            ref_loci=title.split('|')[1]
            dict1[fa_loci]=ref_loci
    return dict1

def merge_list(filepath,suffix,g_folder,inputfolder):
    split=filepath.split(':')
    file_number=len(split)
    mergelist=[]
    for i in range(0,file_number):
        file_path=split[i] + inputfolder +  folder4 + g_folder + '/' + suffix + '/'
        file_list=files_name_listdir_local(file_path,suffix)
        mergelist.append(file_list)
    new_list=list(set(chain.from_iterable(mergelist)))
    return new_list


def concatenate_fa(g_folder,filepath,inputfolder,outputfolder):
    split=filepath.split(':')
    file_number=len(split)
    output_dir = outputfolder + folder5 + g_folder + '/'
    make_dir(output_dir)
    files_list = merge_list(filepath,'final_contig',g_folder,inputfolder)
    for inputfile in files_list:
        outputfile = output_dir + os.path.splitext(inputfile)[0] + '-concatenate.fa'
        fa_dict = fa2dict(filepath,inputfile,'final_contig',g_folder,inputfolder)
        with open(outputfile,'w') as f1:
            for key,value in fa_dict.items():
                f1.writelines("%s\n%s\n" % (key,value))

def concatenate_fa2(g_folder,filepath,inputfolder,outputfolder,ref_file,ref_file2,outgroup):
    split=filepath.split(':')
    file_number=len(split)
    output_dir = outputfolder + folder5 + g_folder + '/'
    make_dir(output_dir)
    files_list = merge_list(filepath,'final_contig',g_folder,inputfolder)
    ref2_dict=fa2dict2(ref_file2)
    with open(ref_file,'r') as f0:
        ref_dict=dict((title,seq) for title,seq in SimpleFastaParser(f0))
    with open(ref_file2,'r') as f2:
        ref_dict2=dict((title,seq) for title,seq in SimpleFastaParser(f2))
    with open(outgroup,'r') as f3:
        ref_dict3=dict((title,seq) for title,seq in SimpleFastaParser(f3))
    for inputfile in files_list:
        outputfile = output_dir + os.path.splitext(inputfile)[0] + '-concatenate.fa'
        fa_dict = fa2dict(filepath,inputfile,'final_contig',g_folder,inputfolder)
        with open(outputfile,'w') as f1:
            title=os.path.splitext(os.path.splitext(inputfile)[0])[0].replace('-','|')
            fa_loci=title.split('|')[0]
            out_seq=ref_dict3[fa_loci]
            f1.writelines(">%s\n%s\n" % (fa_loci,out_seq))
            seq=ref_dict[title]
            f1.writelines(">%s\n%s\n" % (title,seq))
            ref_loci=ref2_dict[fa_loci]
            ref2_title=fa_loci+'|'+ref_loci
            seq2=ref_dict2[ref2_title]
            f1.writelines(">%s\n%s\n" % (ref2_title,seq2))
            for key,value in fa_dict.items():
                f1.writelines("%s\n%s\n" % (key,value))

def all_species(inputfile,species_split):
    split= species_split.split(':')
    fa_list=[]
    fa_dict={}
    num=0
    with open(inputfile,'r') as f1:
        lines=f1.readlines()
        for line in lines:
            if line.startswith('>'):
                name=line.rstrip().split()[-1]
                fa_dict[name]=''
                if name not in fa_list:
                    fa_list.append(name)
            else:
                fa_dict[name]+=line
        for i in split:
            j=i[:-1]
            if j in fa_list:
                num+=1
    return num

def mafft_align(g_folder,outputfolder,species_split,common_num):
    input_dir = outputfolder + folder5 + g_folder +'/'
    output_dir = outputfolder + folder6 + g_folder + '/'
    make_dir(output_dir)
    fa_list = files_name_listdir_local(input_dir, 'fa')
    for fasta in fa_list:
        inputfa = input_dir + os.path.basename(fasta)
        identity = all_species(inputfa,species_split)
        if identity >= int(common_num) :
            outputfa = output_dir + os.path.splitext(fasta)[0]+'.mafft_align'  
            mafft_cline = "mafft " + inputfa + " > " + outputfa
            subp = subprocess.Popen(str(mafft_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
            subp.wait()    
            if subp.poll ==0:
                print("mafft align OK")
            else:
                print("mafft align failed")

def call_consensus(g_folder,outputfolder):
    inputfolder = outputfolder + folder6 + g_folder +'/'
    output_dir = outputfolder + folder7 + g_folder + '/'
    make_dir(output_dir)
    fa_list = files_name_listdir_local(inputfolder, 'mafft_align')
    for fa_file in fa_list:
        alignment = AlignIO.read(inputfolder + fa_file,'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus()
        outputfile= output_dir + os.path.splitext(fa_file)[0] + '.consensus'
        with open(outputfile,'w') as f:
            title= "consensus-" + os.path.splitext(os.path.splitext(fa_file)[0])[0]
            f.writelines(">%s\n%s\n" % (title,consensus))


def call_consensus2(g_folder,outputfolder):
    inputfolder = outputfolder + folder6 + g_folder +'/'
    output_dir = outputfolder + folder7 + g_folder + '/'
    make_dir(output_dir)
    fa_list = files_name_listdir_local(inputfolder, 'mafft_align')
    for fa_file in fa_list:
        fa_dict={}
        outputfile= output_dir + os.path.splitext(fa_file)[0] + '.consensus'
        with open(inputfolder+fa_file,'r') as f1, open(outputfile,'w') as f2:
             for title,seq in SimpleFastaParser(f1):
                 fa_dict[title]=seq.replace('\n','').upper()
             key_list=[key for key in fa_dict.keys()]
             species_num=len(fa_dict.keys())
             seq_len=len(fa_dict[key_list[0]])
             consensus_seq=''
             for i in range(0,seq_len):
                 base_list=[]
                 for j in range(0,species_num):
                     base_list.append(fa_dict[key_list[j]][i])
                 Count_result=Counter(base_list).most_common(2)
                 common_list=[i[0] for i in Count_result if i[0] in 'ATCG' ]
                 consensus_seq+=common_list[0]
             consensus_title= "consensus-" + os.path.splitext(os.path.splitext(fa_file)[0])[0]
             f2.writelines(">%s\n%s\n" % (consensus_title,consensus_seq))   
           

def cp_file(inputfolder,species_split,fq_folder,outputfolder):
    mafft_folder = outputfolder + folder6 + fq_folder +'/'
    file_list= files_name_listdir_local(mafft_folder, 'mafft_align')
    split=species_split.split(':')
    for i in range(0,len(split)):
        species_name=split[i]
        folder_path = species_name + inputfolder + folder3 + fq_folder+'/'
        output_path = outputfolder + folder8 + species_name+'/'+ fq_folder
        make_dir(output_path)
        for i in file_list:
            prefix=i.split('.txt')[0]
            fq1_name=folder_path+prefix+'.txt_R1.fq'
            fq2_name=folder_path+prefix+'.txt_R2.fq'
            if os.path.exists(fq1_name):
                shutil.copy(fq1_name,output_path)
                shutil.copy(fq2_name,output_path)
            else:
                print("%s not found\n%s not found" % (fq1_name,fq2_name))
                pass

def change_title(g_folder,outputfolder,fa_split,ref_number):
    inputfolder = outputfolder + folder6 + g_folder +'/'
    output_dir = outputfolder + folder7 + g_folder + '/'
    make_dir(output_dir)
    fa_list = files_name_listdir_local(inputfolder, 'mafft_align')
    split=fa_split.split(':')
    file_number=len(split)
    for fa_file in fa_list:
        outputfile = output_dir + os.path.splitext(fa_file)[0] + '.fas'
        with open(inputfolder+fa_file,'r') as f1:
            fa_dict={}
            lines=f1.readlines()
            for line in lines:
                num=len(fa_dict.keys())
                if num < int(ref_number):
                    if line.startswith('>'):
                        name=split[num]
                        fa_dict[name]=''    
                    else:
                        fa_dict[name]+=line.replace('\n','').upper()
                else:
                    if line.startswith('>contig1') or line.startswith('>contig2'): 
                        name=line.split(']')[-1].strip()
                        fa_dict[name]=''
                    else:
                        fa_dict[name]+=line.replace('\n','').upper()
            for element in split:
                if element not in fa_dict.keys():
                    seq_len=len(fa_dict[split[0]])
                    fa_dict[element] = '-' * seq_len
            with open(outputfile,'w') as f2:
                for key,value in fa_dict.items():
                    f2.writelines(">%s\n%s\n" % (key,value))
   
def merge_fa(g_folder,outputfolder,fa_split,phy_file,nex_file):
    input_dir = outputfolder + folder7 + g_folder +'/'
    output_dir = outputfolder + folder8 + g_folder + '/'
    make_dir(output_dir)
    fa_list = files_name_listdir_local(input_dir,'fas')
    split=fa_split.split(':')
    file_number=len(split)
    fa_dict={}
    nex_list=[]    
    for fasta in fa_list:
        inputfa = input_dir + fasta
        with open(inputfa,'r') as f0:
            for title,seq in SimpleFastaParser(f0):
                if title not in fa_dict.keys():
                    fa_dict[title]=seq
                else:       
                    fa_dict[title]+=seq                         
    with open(output_dir+phy_file,'w') as f1:
        name=split[1]
        fa_length=len(fa_dict[name])
        f1.writelines("%s\t%s\n" % (file_number,fa_length))       
        for key,value in fa_dict.items():
            f1.writelines("%s\t%s\n" % (key,value))
    with open(output_dir + nex_file,'w') as f2:
        fa_list.sort()
        for fasta in fa_list:
            inputfa = input_dir + fasta
            loci=fasta.split('-')[0]
            with open(inputfa,'r') as f0:
                for title,seq in SimpleFastaParser(f0):
                    if nex_list == []:
                        start=1
                        seq_length=len(seq)
                        end=start + seq_length - 1
                        line='['+loci+' '+str(seq_length)+'bp '+str(start)+'-'+str(end)+' '+'bp]'
                        f2.writelines("%s\n%s\t%s\n" % (line,title,seq))
                        start=end
                        nex_list.append(loci)
                    else: 
                        if loci not in nex_list:
                            seq_length=len(seq)
                            end=start + seq_length - 1
                            line='['+loci+' '+str(seq_length)+'bp '+str(start)+'-'+str(end)+' '+'bp]'
                            f2.writelines("\n%s\n%s\t%s\n" % (line,title,seq))
                            start=end
                            nex_list.append(loci)
                        else:
                            f2.writelines("%s\t%s\n" % (title,seq))
                       
        
if __name__=="__main__":
    concatenate_fa(args.g1folder,args.species_split,args.inputfolder,args.outputfolder)
    concatenate_fa(args.g2folder,args.species_split,args.inputfolder,args.outputfolder)
    #concatenate_fa2(args.g1folder,args.species_split,args.inputfolder,args.outputfolder,args.g1ref,args.g1ref2,args.outgroup)
    #concatenate_fa2(args.g2folder,args.species_split,args.inputfolder,args.outputfolder,args.g2ref,args.g2ref2,args.outgroup)
    mafft_align(args.g1folder,args.outputfolder,args.species_split,args.mafft_number)
    mafft_align(args.g2folder,args.outputfolder,args.species_split,args.mafft_number)
    call_consensus2(args.g1folder,args.outputfolder)
    call_consensus2(args.g2folder,args.outputfolder)
    #cp_file(args.inputfolder,args.species_split,args.g1folder,args.outputfolder)
    #cp_file(args.inputfolder,args.species_split,args.g2folder,args.outputfolder)
    #change_title(args.g1folder,args.outputfolder,args.fa_split,args.ref_number)
    #change_title(args.g2folder,args.outputfolder,args.fa_split,args.ref_number)
    #merge_fa(args.g1folder,args.outputfolder,args.fa_split,args.phy_file1,args.nex_file1)
    #merge_fa(args.g2folder,args.outputfolder,args.fa_split,args.phy_file2,args.nex_file2)
