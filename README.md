# Sequence Alignment and Annotation

Uses Snakemake pipeline for sequence alignment and annotation. Needs Snakemake environments to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Complete Pipeline Instructions

1. Copy sequence files to folder ```data```
2. Rename the sequence files to follow the following naming convention:
```
sample_1.fastq #forward read
sample_2.fastq #reverse read
sample.fastq #long read
``` 
3. Open ```Snakefile``` in a text editor and edit the value for sample. If multiple samples are to be analysed, format it as such: ```[sample1,sample2]```
4. In the ```snakemake``` conda environment, run ```snakemake -c 4 --use-conda```
5. Results will be output in ```results``` folder.