# Mapper_eva
### Code for paper "X-Mapper: fast and accurate sequence alignment via x-mers"

X-Mapper is a fast and accurate short read alignment tool via x-mers.

### Scripts to evaluate alignment performance
- `Alignment_accuracy.ipynb`: Code to evaluate alignment accuracy and run time.
- `Alignment_accuracy_ANI.ipynb`: Code to evaluate alignment accuracy for genomes with various ANI and run time.
- `Alignment_consistency.ipynb`: Code to evaluate alignment consistency.
- `Folder examples`: Scripts and results supporting Fig. S2

**Generating code for evaluating alignment accuracy**
```
python accuracy/SNP_modelpenalty.py
python accuracy/SNP_modelpenalty_ANI.py
```
**Generating code for evaluating alignment accuracy**
```
python consistency/genome_multi_prepare.py
python consistency/SNP_modelmulti_genome.py
python consistency/samfiltersecondary.py
```
**Generating code for evaluating alignment run time**
```
python runtime/SNP_model_time.py
python runtime/timesum_parallel_mapper.py
```
