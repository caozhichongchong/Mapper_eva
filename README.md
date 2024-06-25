# Mapper_eva
### Code for paper "Mapper: fast and accurate sequence alignment via x-mers"

Mapper is a fast and accurate short read alignment tool via x-mers.

### Scripts to evaluate alignment performance
- `Alignment_accuracy.ipynb`: Code to evaluate alignment accuracy and run time.
- `Alignment_consistency.ipynb`: Code to evaluate alignment consistency.
- `Folder examples`: Scripts and results supporting Figs. S1-S3

**Generating code for evaluating alignment accuracy**
```
python accuracy/SNP_modelpenalty.py
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
