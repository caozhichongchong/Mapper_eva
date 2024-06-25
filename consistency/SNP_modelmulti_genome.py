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
                      type=str, default='/home/caozhichongchong/WS/WS1/genome/donor_species/SNP_curate/test_data/multi_Genome/',
                      metavar='input/')
# optional output setup
optional.add_argument("-fa",
                      help="file extension of fasta files",
                      type=str, default='_final.scaffolds.fasta',
                      metavar='.corrected.fasta')
optional.add_argument("-fq",
                      help="file extension of fastq files",
                      type=str, default='_1.fastq',
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
include_secondary = True
prep = False
# penalty: snp, indel start, indel extension, match, ambiguity
penaltyset = {
    'bowtie':[6,5,3,0,1],
    'minimap':[6,4,2,0,1],
    'bwa':[5,6,1,0,1],
    'mapper':[1*10,2*10,int(0.5*10),0*10,int(0.1*10)]
}
input_script = args.s
genome_root = args.i
output_dir = args.o + '/SNP_model_multi'
genome_name = args.fa
fastq_name=args.fq
input_script_sub = '%s/SNP_model_multi'%(input_script)
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
    cmds = ''
    if 'mapper1' not in sam:
        cmds = 'samtools view -F 4 %s > %s.mapped\n'%(sam,sam)# remove unmapped alignments
        cmds += 'mv %s.mapped %s\n'%(sam,sam)
    cmds += 'conda deactivate\npython /home/caozhichongchong/WS/WS2/github/snp_finder/snp_finder/scripts/samfiltersecondary.py -i %s\n'%(sam)
    if 'mapper1' not in sam:
        cmds += 'conda activate bt\n'
    return cmds

def run_bowtie(files,database,tempbamoutput):
    #
    # generate code
    cmds = 'bowtie2-build %s %s\n'%(database,database)
    #-a mode: search for and report all alignments, including secondary
    #-k mode: search for one or more alignments, report each
    # not using --ma for end to end model, --ignore-quals for mismatch quality, always use the highest penalty
    cmds += '/usr/bin/time -v bowtie2 --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals -x %s -1 %s -2 %s -S %s.sam\n' % (
            min(40, args.t), penalty[0],penalty[0],penalty[1],penalty[2],penalty[1],penalty[2],penalty[4], database, files, files.replace('_1.fastq','_2.fastq'),tempbamoutput)
    if check_bowtieminimap_kmer:
        cmds += '/usr/bin/time -v bowtie2 --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals --seedlen=15 -x %s -1 %s -2 %s -S %s.15mer.sam\n' % (
            min(40, args.t), penalty[0], penalty[0], penalty[1], penalty[2], penalty[1], penalty[2], penalty[4], database,
            files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.15mer.sam')
    if include_secondary:
        cmds += '/usr/bin/time -v bowtie2 -a --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals -x %s -1 %s -2 %s -S %s.all.sam\n' % (
            min(40, args.t), penalty[0], penalty[0], penalty[1], penalty[2], penalty[1], penalty[2], penalty[4],
            database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.all.sam')
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    return cmds

def run_minimap(files,database,tempbamoutput):
    # default setting claims to retain all primary mappings. -N output top secondary mappings, if higher than -p
    cmds = 'minimap2 -d %s.mmi %s \n' % (database, database)
    # -A match score set as 2, -B -= 2 and --score-N += 2
    if penalty[0] > 2:
        cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s %s >%s.sam\n' % (
                min(40, args.t), penalty[0] - 2,penalty[1] - 2,penalty[2],penalty[3] + 2 ,penalty[3]-penalty[4] + 2, database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
    else:
        cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s %s >%s.sam\n' % (
                min(40, args.t), penalty[0] - 0.5,penalty[1] - 0.5,penalty[2],penalty[3] + 0.5 ,penalty[3]-penalty[4] + 0.5, database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    if check_bowtieminimap_kmer:
        # -A match score set as 2, -B -= 2 and --score-N += 2
        if penalty[0] > 2:
            cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s -k10 %s.mmi %s %s >%s.10mer.sam\n' % (
                min(40, args.t), penalty[0] - 2, penalty[1] - 2, penalty[2], penalty[3] + 2, penalty[3] - penalty[4] + 2,
                database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
        else:
            cmds += '/usr/bin/time -v minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s -k10 %s.mmi %s %s >%s.10mer.sam\n' % (
                min(40, args.t), penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5, penalty[3] - penalty[4] + 0.5,
                database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.10mer.sam')
    if include_secondary:
        if penalty[0] > 2:
            cmds += '/usr/bin/time -v minimap2 -ax sr -N 100 -t %s -B %s -O %s -E %s  -A %s --score-N %s %s.mmi %s %s >%s.all.sam\n' % (
                min(40, args.t), penalty[0] - 2, penalty[1] - 2, penalty[2], penalty[3] + 2, penalty[3] - penalty[4] + 2,
                database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
        else:
            cmds += '/usr/bin/time -v minimap2 -ax sr -N 100 -t %s -B %s -O %s -E %s  -A %s --score-N %s %s.mmi %s %s >%s.all.sam\n' % (
                min(40, args.t), penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5,
                penalty[3] - penalty[4] + 0.5, database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.all.sam')
    return cmds

def run_bwa(files,database,tempbamoutput):
    # -a output all alignments for SE or unpaired PE
    cmds = 'bwa index %s\n'%(database)
    # -A match score set as 2, -B -= 2
    if penalty[0] > 2:
        cmds += '/usr/bin/time -v bwa mem -t %s -B %s -O %s -E %s -A %s %s %s %s > %s.sam\n' % (
                min(40, args.t), penalty[0] - 2,penalty[1] - 2,penalty[2],penalty[3]+ 2, database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
    else:
        cmds += '/usr/bin/time -v bwa mem -t %s -B %s -O %s -E %s -A %s %s %s %s > %s.sam\n' % (
            min(40, args.t), penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5, database, files, files.replace('_1.fastq','_2.fastq'),
            tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    if include_secondary:
        if penalty[0] > 2:
            cmds += '/usr/bin/time -v bwa mem -a -t %s -B %s -O %s -E %s -A %s %s %s %s > %s.all.sam\n' % (
                min(40, args.t), penalty[0] - 2, penalty[1] - 2, penalty[2], penalty[3] + 2, database, files, files.replace('_1.fastq','_2.fastq'), tempbamoutput)
        else:
            cmds += '/usr/bin/time -v bwa mem -a -t %s -B %s -O %s -E %s -A %s %s %s %s > %s.all.sam\n' % (
                min(40, args.t), penalty[0] - 0.5, penalty[1] - 0.5, penalty[2], penalty[3] + 0.5, database, files, files.replace('_1.fastq','_2.fastq'),
                tempbamoutput)
        cmds += remove_unmappedreads(tempbamoutput + '.all.sam')
    return cmds

def run_mapper(files,database,tempbamoutput):
    penalty2 = [penalty[3]+penalty[0],penalty[1],penalty[2],penalty[-1]]
    max_penalty = penalty2[0]*0.1
    cmds = ''
    # spacing penalty <expected> <distancePerPenalty>  distancePerPenalty default = 50, new penalty = default/SNP_penalty
    cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --max-penalty-span 0 --no-infer-ancestors --num-threads %s --reference %s --paired-queries %s %s --spacing 100 %s  --out-sam %s.sam\n' % (
                latest_mapper, max_penalty, penalty2[0],penalty2[1],penalty2[2],penalty2[3],min(40, args.t),database, files, files.replace('_1.fastq','_2.fastq'), int(50/penalty2[0]), tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    # cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --max-penalty-span 0 --num-threads %s --reference %s --paired-queries %s %s --spacing 100 %s  --out-sam %s.ancestor.sam\n' % (
    #     latest_mapper, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database,
    #     files, files.replace('_1.fastq', '_2.fastq'), int(50 / penalty2[0]), tempbamoutput)
    # cmds += remove_unmappedreads(tempbamoutput + '.ancestor.sam')
    # cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --no-gapmers --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --max-penalty-span 0 --no-infer-ancestors --num-threads %s --reference %s --paired-queries %s %s --spacing 100 %s  --out-sam %s.nogap.sam\n' % (
    #     latest_mapper, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database,
    #     files, files.replace('_1.fastq', '_2.fastq'), int(50 / penalty2[0]), tempbamoutput)
    # cmds += remove_unmappedreads(tempbamoutput + '.nogap.sam')
    # cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --no-infer-ancestors  --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --num-threads %s --reference %s --paired-queries %s %s   --out-sam %s.noancestor.sam\n' % (
    #     latest_mapper, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database, files, files.replace('_1.fastq','_2.fastq'),
    #     tempbamoutput)
    # cmds += remove_unmappedreads(tempbamoutput + '.noancestor.sam')
    # kmer default: minimap 21, bowtie 20-22 for --sensitive mode, bwa 19
    cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --block-length 12 --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --max-penalty-span 0 --no-infer-ancestors --num-threads %s --reference %s --paired-queries %s %s --spacing 100 %s  --out-sam %s.kmer12.sam\n' % (
        latest_mapper_kmer, max_penalty, penalty2[0],penalty2[1],penalty2[2],penalty2[3],min(40, args.t),database, files, files.replace('_1.fastq','_2.fastq'), int(50/penalty2[0]), tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.kmer12.sam')
    cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --block-length 18 --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --max-penalty-span 0 --no-infer-ancestors --num-threads %s --reference %s --paired-queries %s %s --spacing 100 %s  --out-sam %s.kmer18.sam\n' % (
        latest_mapper_kmer, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database,
        files, files.replace('_1.fastq', '_2.fastq'), int(50 / penalty2[0]), tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.kmer18.sam')
    cmds += '/usr/bin/time -v java -Xms55g -Xmx55g -jar %s --block-length 24 --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --ambiguity-penalty %s --max-penalty-span 0 --no-infer-ancestors --num-threads %s --reference %s --paired-queries %s %s --spacing 100 %s  --out-sam %s.kmer24.sam\n' % (
        latest_mapper_kmer, max_penalty, penalty2[0], penalty2[1], penalty2[2], penalty2[3], min(40, args.t), database,
        files, files.replace('_1.fastq', '_2.fastq'), int(50 / penalty2[0]), tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.kmer24.sam')
    if tool == 'bwa':
        cmds = cmds.replace('--ambiguity-penalty %s '%(penalty2[3]),'')
    return cmds

################################################## Main ########################################################
# find all files
allfolders = glob.glob('%s/*'%(genome_root))
for folder in allfolders:
    allgenome = glob.glob('%s/*%s'%(folder,fastq_name))
    if prep:
        for fastq_file in allgenome:
            os.system('#cat %s %s > %s'%(fastq_file,fastq_file.replace('_1.fastq','_2.fastq'),
                                        fastq_file.replace('_1.fastq','_all.fastq')))
        os.system('#cat %s/*_all.fastq > %s/all.fastq'%(folder,folder))
        os.system('cat %s/*_final.scaffolds.fasta > %s/final.scaffolds.fasta' % (folder, folder))
    fastq_name = '_1.fastq'
    allgenome = glob.glob('%s/*%s'%(folder,fastq_name))
    print(allgenome)
    for fastq_file in allgenome:
        databaseset = ['%s/%s'%(folder,os.path.split(fastq_file)[-1].replace(fastq_name,genome_name))]
        print(fastq_file)
        for database in databaseset:
            print(database)
            sample_name = os.path.split(fastq_file)[-1]
            database_name = os.path.split(folder)[-1]
            # find fastq
            for tool in penaltyset:
                penalty = penaltyset[tool]
                # mapper
                tempbamoutput = '%s/%s/%s_%s.mapper1' % (output_dir, tool, sample_name,database_name)
                cmds = '#!/bin/bash\nsource ~/.bashrc\n'
                cmds += run_mapper(fastq_file, database, tempbamoutput)
                # # minimap
                # cmds += 'conda activate bt\n'
                # tempbamoutput = '%s/%s/%s_%s.minimap' % (output_dir, tool, sample_name,database_name)
                # cmds += run_minimap(fastq_file, database, tempbamoutput)
                # # bwa
                # tempbamoutput = '%s/%s/%s_%s.bwa' % (output_dir, tool, sample_name,database_name)
                # cmds += run_bwa(fastq_file, database, tempbamoutput)
                # # bowtie
                # tempbamoutput = '%s/%s/%s_%s.bowtie'%(output_dir,tool,sample_name,database_name)
                # cmds += run_bowtie(fastq_file, database, tempbamoutput)
                f1 = open(os.path.join(input_script_sub, '%s.%s.%s.sh' % (sample_name,database_name,tool)), 'w')
                f1.write(cmds)
                f1.close()

f1 = open(os.path.join(input_script, 'allmultitest.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
    f1.write('bash %s 2> %s.err 1> %s.out \n' % (sub_scripts, sub_scripts,sub_scripts))
f1.close()

################################################### END ########################################################
