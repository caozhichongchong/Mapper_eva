# step 1 modelling SNPs and test 2 methods
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS of each species",
                      type=str, default='/home/caozhichongchong/WS/WS1/genome/donor_species/SNP_curate/test_data/',
                      metavar='input/')
# optional output setup
optional.add_argument("-fa",
                      help="file extension of fasta files",
                      type=str, default='.fasta.corrected',
                      metavar='.fasta')
optional.add_argument("-fq",
                      help="file extension of fastq files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')

optional.add_argument("-s",
                      help="a folder for your mapper",
                      type=str, default='/home/caozhichongchong/WS/WS2/scripts/snp_curate/',
                      metavar='mapper_folder/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='/home/caozhichongchong/WS/WS1/genome/donor_species/SNP_curate/SNP_model_parallel/',
                      metavar='.')
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 30)",
                      metavar="1 or more", action='store', default=30, type=int)
optional.add_argument("-indel",
                      help="whether to insert indels",
                      type=str, default='True',
                      metavar='False or True')

################################################## Definition ########################################################
args = parser.parse_args()
latest_mapper = glob.glob('%s/mapper-1*.jar'%(args.s))
latest_mapper.sort()
latest_mapper=latest_mapper[-1]
threadstouse=args.t
output_dir = args.o
input_script = args.s
try:
    os.mkdir(output_dir)
except IOError:
    pass
input_script_sub = '%s/SNP_model_parallel_new'%(input_script)
try:
    os.mkdir(input_script_sub)
except IOError:
    pass
def run_bowtie(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'conda activate bt\n'
    for i in range(1, 30, 2):
        if 'all_final' not in database:
            cmds += 'time timeout 1h bash -c "bowtie2-build %s %s && bowtie2 --threads %s -x %s -1 %s -2 %s -S %s.sam"\n' % (
                database, database, i, database, files, files2, tempbamoutput)
        else:
            #all alignment
            cmds += 'time timeout 1h bash -c "bowtie2-build %s %s && bowtie2 -a --threads %s -x %s -1 %s -2 %s -S %s.sam"\n' % (
                database, database, i, database, files, files2, tempbamoutput)
    return cmds

def run_minimap(files,files2,database,tempbamoutput):
    cmds = 'conda activate bt\n'
    for i in range(1, 30, 2):
        if 'all_final' not in database:
            cmds += 'time timeout 1h bash -c "minimap2 -t %s -d %s.mmi %s && minimap2 -ax sr -t %s %s.mmi %s %s >%s.sam"\n' % (
                i,
                database, database, i,
                database, files, files2, tempbamoutput)
        else:
            cmds += 'time timeout 1h bash -c "minimap2 -t %s -d %s.mmi %s && minimap2 -ax sr -N 100 -t %s %s.mmi %s %s >%s.sam"\n' % (
                i,
                database, database, i,
                database, files, files2, tempbamoutput)
    return cmds

def run_mapper(files,files2,database,tempbamoutput):
    cmds = ''
    for i in range(1, 30, 2):
        cmds += 'time timeout 1h bash -c "java -Xms30g -Xmx30g -jar %s --no-infer-ancestors --num-threads %s --reference %s --paired-queries %s %s --spacing 100 50  --out-sam %s.sam"\n' % (
            latest_mapper, i, database, files, files2, tempbamoutput)
    return cmds

def run_bwa(files,files2,database,tempbamoutput):
    cmds = 'conda activate bt\n'
    for i in range(1, 30, 2):
        if 'all_final' not in database:
            cmds += 'time timeout 1h bash -c "bwa index %s && bwa mem -t %s %s %s %s > %s.sam"\n' % (
                database, i, database, files, files2, tempbamoutput)
        else:
            cmds += 'time timeout 1h bash -c "bwa index %s && bwa mem -a -t %s %s %s %s > %s.sam"\n' % (
                database, i, database, files, files2, tempbamoutput)
    return cmds

databaseset = glob.glob('%s/new/am_BaFr_g0050.fasta'%(args.i)) + \
glob.glob('%s/penalty_test/am0230.fasta'%(args.i))+ \
glob.glob('%s/multi_Genome/Bactor/all_final.scaffolds.fasta'%(args.i))
print(databaseset)
for database in databaseset:
    file1 = database.replace('_final.scaffolds.fasta','_1.fastq').replace('.fasta','_1.fastq')
    file2 = file1.replace('_1.fastq','_2.fastq')
    print(database,file1,file2)
    databasename = os.path.basename(file1).split('_1.fastq')[0]
    # call SNPs by time bowtie2
    cmds = ''
    cmds += run_bowtie(file1, file2,
                       database,
                       os.path.join(output_dir,
                                    databasename + '.bowtie'))
    f1 = open(os.path.join(input_script_sub, '%s.bowtie.sh' % (databasename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()
    # call SNPs by time bwa
    cmds = ''
    cmds += run_bwa(file1, file2,
                    database,
                    os.path.join(output_dir,
                                 databasename + '.bwa'))
    f1 = open(os.path.join(input_script_sub, '%s.bwa.sh' % (databasename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()
    # call SNPs by time minimap
    cmds = ''
    cmds += run_minimap(file1, file2,
                        database,
                        os.path.join(output_dir,
                                     databasename + '.minimap'))
    f1 = open(os.path.join(input_script_sub, '%s.minimap.sh' % (databasename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()
    # call SNPs by time mapper
    cmds = ''
    cmds += run_mapper(file1, file2,
                       database,
                       os.path.join(output_dir,
                                    databasename + '.mapper'))
    f1 = open(os.path.join(input_script_sub, '%s.mapper.sh' % (databasename)),
              'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()
    # run bowtie and mapper on the same node
    cmds = ''
    cmds += 'bash %s 2> %s.out 1> %s.err\n' % (
        os.path.join(input_script_sub, '%s.mapper.sh' % (databasename)),
        os.path.join(input_script_sub, '%s.mapper.sh' % (databasename)),
        os.path.join(input_script_sub, '%s.mapper.sh' % (databasename)))
    cmds += 'bash %s 2> %s.out 1> %s.err\n' % (
        os.path.join(input_script_sub, '%s.bowtie.sh' % (databasename)),
    os.path.join(input_script_sub, '%s.bowtie.sh' % (databasename)),
    os.path.join(input_script_sub, '%s.bowtie.sh' % (databasename)))
    cmds += 'bash %s 2> %s.out 1> %s.err\n' % (
        os.path.join(input_script_sub, '%s.minimap.sh' % (databasename)),
        os.path.join(input_script_sub, '%s.minimap.sh' % (databasename)),
        os.path.join(input_script_sub, '%s.minimap.sh' % (databasename))
    )
    cmds += 'bash %s 2> %s.out 1> %s.err\n' % (
        os.path.join(input_script_sub, '%s.bwa.sh' % (databasename)),
    os.path.join(input_script_sub, '%s.bwa.sh' % (databasename)),
    os.path.join(input_script_sub, '%s.bwa.sh' % (databasename)))
    f1 = open(os.path.join(input_script_sub, '%s.all.vcf.sh0' % (databasename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()