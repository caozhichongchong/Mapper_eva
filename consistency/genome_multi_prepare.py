import glob,os
from Bio import SeqIO
from Bio.Seq import Seq

inputfolder = '/home/caozhichongchong/WS/WS1/genome/donor_species/SNP_curate/test_data/multi_Genome/Bactor/allfasta'
allfasta = glob.glob('%s/*fasta'%(inputfolder))

allname = []
for fasta in allfasta:
    fasta_name = '_'.join(os.path.basename(fasta).split('_')[:2])
    print(fasta_name)
    if fasta_name in allname:
        print('duplicatename',fasta_name)
    else:
        allname.append(fasta_name)
    newseq = []
    for record in SeqIO.parse(fasta,'fasta'):
        if fasta_name in str(record.id):
            break
        newseq.append('>%s_%s\n%s\n'%(fasta_name,str(record.id),str(record.seq)))
    if newseq != []:
        f1 = open(fasta + '.rename','w')
        f1.write(''.join(newseq))
        f1.close()

inputfolder = '/home/caozhichongchong/WS/WS1/genome/donor_species/SNP_curate/test_data/multi_Genome/Bifido/allfasta'
allfasta = glob.glob('%s/*fasta'%(inputfolder))

allname = []
for fasta in allfasta:
    fasta_name = '_'.join(os.path.basename(fasta).split('_')[:2])
    print(fasta_name)
    if fasta_name in allname:
        print('duplicatename',fasta_name)
    else:
        allname.append(fasta_name)
    newseq = []
    for record in SeqIO.parse(fasta,'fasta'):
        if fasta_name in str(record.id):
            break
        newseq.append('>%s_%s\n%s\n'%(fasta_name,str(record.id),str(record.seq)))
    if newseq != []:
        f1 = open(fasta + '.rename','w')
        f1.write(''.join(newseq))
        f1.close()

