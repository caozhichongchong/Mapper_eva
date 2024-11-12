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
                      type=str, default='/media/caozhichongchong/0fca1a13-dd3f-43f5-9b94-876be77c58da/WS1/genome/donor_species/SNP_curate/test_data/penalty_test/',
                      metavar='input/')
# optional output setup
optional.add_argument("-fa",
                      help="file extension of fasta files",
                      type=str, default='.fasta',
                      metavar='.corrected.fasta')
optional.add_argument("-fq",
                      help="file extension of fastq files",
                      type=str, default='_all.fastq',
                      metavar='_all.fastq')

optional.add_argument("-s",
                      help="a folder for your mapper",
                      type=str, default='/home/caozhichongchong/WS/WS2/scripts/snp_curate/',
                      metavar='mapper_folder/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='/media/caozhichongchong/0fca1a13-dd3f-43f5-9b94-876be77c58da/WS1/genome/donor_species/SNP_curate/',
                      metavar='.')
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 30)",
                      metavar="1 or more", action='store', default=30, type=int)

################################################## Definition ########################################################
# setting match score as 0 for all
args = parser.parse_args()
check_bowtieminimap_kmer = False
include_secondary = False
# penalty: snp, indel start, indel extension, match, ambiguity
penaltyset = {
    'bowtie':[6,5,3,0,1],
    'minimap':[6,4,2,0,1],
    'bwa':[5,6,1,0,1],
    'mapper':[1*10,2*10,int(0.5*10),0*10,int(0.1*10)]
}
input_script = args.s
genome_root = args.i
output_dir = args.o + '/SNP_model_penalty'
genome_name = args.fa
fastq_name=args.fq
input_script_sub = '%s/SNP_model_penalty'%(input_script)
latest_mapper = glob.glob('%s/mapper-1*.jar'%(args.s))
latest_mapper = [x for x in latest_mapper if 'experimental' not in x][0]
print(latest_mapper)
latest_mapper_kmer = glob.glob('%s/mapper-1*kmers*.jar'%(args.s))
latest_mapper_kmer.sort()
latest_mapper_kmer=latest_mapper_kmer[-1]
print(latest_mapper_kmer)
try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/minimap')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/bowtie')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/mapper')
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

def remove_unmappedreads(sam):
    cmds = 'samtools view -F 4,256 %s > %s.mapped\n'%(sam,sam)# remove unmapped and secondary alignments
    cmds += 'mv %s.mapped %s\n'%(sam,sam)
    return cmds

def run_bowtie(files,database,tempbamoutput):
    # generate code
    cmds = 'conda activate bt\n#bowtie2-build %s %s\n'%(database,database)
    # Segmentation fault (core dumped) for bowtie2 -a for Bactor_all
    # not using --ma for end to end model, --ignore-quals for mismatch quality, always use the highest penalty
    cmds += '/usr/bin/time -v bowtie2 --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s -N 1 --ignore-quals -x %s -U %s -S %sN1.sam\n' % (
            min(40, args.t), penalty[0],penalty[0],penalty[1],penalty[2],penalty[1],penalty[2],penalty[4], database, files, tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + 'N1.sam')
    # cmds += '/usr/bin/time -v bowtie2 --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals -x %s -U %s -S %s.sam\n' % (
    #     min(40, args.t), penalty[0], penalty[0], penalty[1], penalty[2], penalty[1], penalty[2], penalty[4], database,
    #     files, tempbamoutput)
    # cmds += remove_unmappedreads(tempbamoutput + '.sam')
    # if check_bowtieminimap_kmer:
    #     cmds += '/usr/bin/time -v bowtie2 --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals --seedlen=12 -x %s -U %s -S %s.12mer.sam\n' % (
    #         min(40, args.t), penalty[0], penalty[0], penalty[1], penalty[2], penalty[1], penalty[2], penalty[4], database,
    #         files, tempbamoutput)
    #     cmds += remove_unmappedreads(tempbamoutput + '.12mer.sam')
    # if include_secondary:
    #     cmds += '/usr/bin/time -v bowtie2 -a --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals -x %s -U %s -S %s.all.sam\n' % (
    #         min(40, args.t), penalty[0], penalty[0], penalty[1], penalty[2], penalty[1], penalty[2], penalty[4],
    #         database, files, tempbamoutput)
    #     cmds += remove_unmappedreads(tempbamoutput + '.all.sam')
    return cmds

def run_minimap(files,database,tempbamoutput):
    cmds = 'minimap2 -d %s.mmi %s \n' % (database, database)
    # -A match score set as 2, -B -= 2 and --score-N += 2
    if penalty[0] > 2:
        cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s >%s.sam\n' % (
                min(40, args.t), penalty[0] - 2,penalty[1] - 2,penalty[2],penalty[3] + 2 ,penalty[3]-penalty[4] + 2, database, files, tempbamoutput)
    else:
        cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s >%s.sam\n' % (
                min(40, args.t), penalty[0] - 0.5,penalty[1] - 0.5,penalty[2] ,penalty[3] + 0.5 ,penalty[3]-penalty[4] + 0.5, database, files, tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    if check_bowtieminimap_kmer:
        cmds += 'minimap2 -d %s.mmi -k10 %s \n' % (database, database)
        # -A match score set as 2, -B -= 2 and --score-N += 2
        if penalty[0] > 2:
            cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s -k12 %s.mmi %s >%s.12mer.sam\n' % (
                min(40, args.t), penalty[0] - 2, penalty[1]- 2, penalty[2], penalty[3] + 2, penalty[3] - penalty[4] + 2,
                database, files, tempbamoutput)
        else:
            cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s -k12 %s.mmi %s >%s.12mer.sam\n' % (
                min(40, args.t), penalty[0] - 0.5, penalty[1]- 0.5, penalty[2], penalty[3] + 0.5, penalty[3] - penalty[4] + 0.5,
                database, files, tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.12mer.sam')
    if include_secondary:
        if penalty[0] > 2:
            cmds += '/usr/bin/time -v minimap2 -ax sr -N 100 -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s >%s.all.sam\n' % (
                min(40, args.t), penalty[0] - 2, penalty[1]- 2, penalty[2], penalty[3] + 2, penalty[3] - penalty[4] + 2,
                database, files, tempbamoutput)
        else:
            cmds += '/usr/bin/time -v minimap2 -ax sr -N 100 -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s >%s.all.sam\n' % (
                min(40, args.t), penalty[0] - 0.5, penalty[1]- 0.5, penalty[2], penalty[3] + 0.5,
                penalty[3] - penalty[4] + 0.5, database, files, tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.all.sam')
    return cmds

def run_bwa(files,database,tempbamoutput):
    cmds = '#bwa index %s\n'%(database)
    # -A match score set as 1, -B -= 1
    if penalty[0] > 2:
        cmds += '/usr/bin/time -v bwa mem -t %s -B %s -O %s -E %s -L 100 -A %s %s %s > %s.sam\n' % (
                min(40, args.t), penalty[0] - 1,penalty[1]-1,penalty[2],penalty[3]+ 2, database, files, tempbamoutput)
    else:
        cmds += '/usr/bin/time -v bwa mem -t %s -B %s -O %s -E %s -L 100 -A %s %s %s > %s.sam\n' % (
            min(40, args.t), penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5, database, files,
            tempbamoutput)
    if include_secondary:
        if penalty[0] > 2:
            cmds += '/usr/bin/time -v bwa mem -a -t %s -B %s -O %s -E %s -L 100 -A %s %s %s > %s.all.sam\n' % (
                min(40, args.t), penalty[0] - 1, penalty[1] -1, penalty[2],  penalty[3] + 2, database, files,  tempbamoutput)
        else:
            cmds += '/usr/bin/time -v bwa mem -a -t %s -B %s -O %s -E %s -L 100 -A %s %s %s > %s.all.sam\n' % (
                min(40, args.t), penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5, database, files,
                tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.all.sam')
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    return cmds

def run_strobealign(files,database,tempbamoutput):
    cmds = 'conda activate strobealign\n'
    cmds += '/usr/bin/time -v strobealign -t %s -B %s -O %s -E %s -A %s -L 1000 %s %s > %s.sam\n'%(
            min(40, args.t),penalty[0],penalty[1],penalty[2],penalty[3],
                                                                                               database,files,tempbamoutput)

    cmds += 'conda activate bt\n'
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    return cmds

def run_last(files,database,tempbamoutput):
    cmds = '#lastdb -P30 -uNEAR %s.mydb %s\n'%(database,database)
    if penalty[0] > 2:
        cmds += '/usr/bin/time -v last-train -P30 -q %s -a %s -b %s -r %s -X 1 -Q1 %s.mydb %s > reads.train\n' % (penalty[0] - 2, penalty[1]- 2, penalty[2], penalty[3] + 2,
                                                                                                                  database, files)
        cmds += '/usr/bin/time -v lastal -P%s -T1 -q %s -a %s -b %s -r %s -X 1 -p reads.train %s.mydb %s > %s.maf\n'%(
                    min(40, args.t),penalty[0] - 2, penalty[1]- 2, penalty[2], penalty[3] + 2,
                                                                                                       database,files,tempbamoutput)
    else:
        cmds += '/usr/bin/time -v last-train -P30 -q %s -a %s -b %s -r %s -X 1 -Q1 %s.mydb %s > reads.train\n' % (
        penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5,
        database, files)
        cmds += '/usr/bin/time -v lastal -P%s -T1 -q %s -a %s -b %s -r %s -X 1 -p reads.train %s.mydb %s > %s.maf\n' % (
            min(40, args.t), penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5,
            database, files, tempbamoutput)
    if 'DRR' in database:
        cmds = cmds.replace(' -X 1 -p ',' -X 1 -K5 -p ')
    cmds += 'maf-convert sam %s.maf > %s.sam\n'%(tempbamoutput,tempbamoutput)
    return cmds

def run_mapper(files,database,tempbamoutput):
    penalty2 = [penalty[3]+penalty[0],penalty[1],penalty[2],penalty[-1]]
    max_penalty = penalty2[0]*0.1
    cmds = ''
    # kmer default: minimap 21, bowtie 20-22 for --sensitive mode, bwa 19
    # # ancestor inference having problems - turned off
    # cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --additional-extend-insertion-penalty 0  --max-penalty-span 0 --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --num-threads %s --reference %s --queries %s  --out-sam %s.withancestor.sam\n' % (
    #             latest_mapper, max_penalty, penalty2[0],penalty2[1],penalty2[2],penalty2[3],min(40, args.t),database, files, tempbamoutput)
    # cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --no-infer-ancestors --additional-extend-insertion-penalty 0  --max-penalty-span 0 --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --num-threads %s --reference %s --queries %s  --out-sam %s.sam\n' % (
    #     latest_mapper, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database, files,
    #     tempbamoutput)
    # cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --no-infer-ancestors --no-gapmers --additional-extend-insertion-penalty 0  --max-penalty-span 0 --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --num-threads %s --reference %s --queries %s  --out-sam %s.nogapped.sam\n' % (
    #     latest_mapper, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database,
    #     files,
    #     tempbamoutput)
    if 'DRR' in database:
        for k in [16]:
            cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --no-gapmers --no-infer-ancestors --additional-extend-insertion-penalty 0  --max-penalty-span 0 --block-length %s --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --num-threads %s --reference %s --queries %s  --out-sam %s.kmer%s.sam\n' % (
                latest_mapper_kmer, k, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database,
                files,
                tempbamoutput,k)
    else:
        for k in [24]:
            cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --no-gapmers --no-infer-ancestors --additional-extend-insertion-penalty 0  --max-penalty-span 0 --block-length %s --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --num-threads %s --reference %s --queries %s  --out-sam %s.kmer%s.sam\n' % (
                latest_mapper_kmer, k, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database,
                files,
                tempbamoutput,k)
    if tool == 'bwa':
        cmds = cmds.replace('--ambiguity-penalty %s '%(penalty2[3]),'')
    return cmds

################################################## Main ########################################################
# find all files
allgenome = glob.glob('%s/*%s'%(genome_root,fastq_name))
print(allgenome)
for fastq_file in allgenome:
    databaseset = ['%s/%s'%(genome_root,os.path.split(fastq_file)[-1].replace(fastq_name,genome_name))]
    for database in databaseset:
        sample_name = os.path.split(fastq_file)[-1]
        database_name = os.path.split(database)[-1]
        # find fastq
        for tool in penaltyset:
            penalty = penaltyset[tool]
            cmds = '#!/bin/bash\nsource ~/.bashrc\n'
            # # mapper
            # tempbamoutput = '%s/%s/%s_%s.mapper1' % (output_dir, tool, sample_name,database_name)
            # cmds += run_mapper(fastq_file, database, tempbamoutput)
            # # minimap
            cmds += 'conda activate bt\n'
            # tempbamoutput = '%s/%s/%s_%s.minimap' % (output_dir, tool, sample_name,database_name)
            # cmds += run_minimap(fastq_file, database, tempbamoutput)
            # bwa
            tempbamoutput = '%s/%s/%s_%s.bwa' % (output_dir, tool, sample_name,database_name)
            cmds += run_bwa(fastq_file, database, tempbamoutput)
            # # bowtie
            # tempbamoutput = '%s/%s/%s_%s.bowtie'%(output_dir,tool,sample_name,database_name)
            # cmds += run_bowtie(fastq_file, database, tempbamoutput)
            # # strobealign
            # tempbamoutput = '%s/%s/%s_%s.strobealign' % (output_dir, tool, sample_name, database_name)
            # cmds += run_strobealign(fastq_file, database, tempbamoutput)
            # # last
            # tempbamoutput = '%s/%s/%s_%s.last' % (output_dir, tool, sample_name, database_name)
            # cmds += run_last(fastq_file, database, tempbamoutput)
            f1 = open(os.path.join(input_script_sub, '%s.%s.%s.sh' % (sample_name,database_name,tool)), 'w')
            f1.write(cmds)
            f1.close()

f1 = open(os.path.join(input_script, 'allpenaltytest.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
    f1.write('bash %s 2> %s.err 1> %s.out \n' % (sub_scripts, sub_scripts,sub_scripts))
f1.close()

################################################### END ########################################################
