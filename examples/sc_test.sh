bwa index  am0230.fasta
minimap2 -ax sr -t 30 -B 4 -O 3 -E 3 -A 2 --score-N 1 am0230.fasta.mmi trial_sc_all.fastq > trial_sc.minimap.sam 
bwa mem -t 30 -B 4 -O 3 -E 3 -A 2 am0230.fasta trial_sc_all.fastq > trial_sc.bwa.sam 
minimap2 -ax sr -t 30 -B 4 -O 3 -E 3 -A 2 --score-N 1 -z 10000,200 am0230.fasta.mmi trial_sc_all.fastq > trial_sc.minimap_nosc.sam
bwa mem -t 30 -B 4 -O 3 -E 3 -A 2 -L 100 am0230.fasta trial_sc_all.fastq > trial_sc.bwa_nosc.sam 
