# Sequence Alignment and Annotation

Uses Snakemake pipeline for sequence alignment and annotation. Needs Snakemake environments to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Complete Pipeline Instructions

1. Create folder ```data``` and copy sequence files to it.
2. Rename the sequence files to follow the following naming convention:
```
sample_1.fastq #forward read
sample_2.fastq #reverse read
sample.fastq #long read
``` 
3. Open ```Snakefile``` in a text editor and edit the value for sample. If multiple samples are to be analysed, format it as such: ```['sample1','sample2']```
4. Depending on whether you want to assemble short reads, long reads or do a hybrid assembly, comment out the ```rule unicycler``` rules you do not need.
5. In the ```snakemake``` conda environment, run ```snakemake --use-conda```.
6. Results will be output in ```results``` folder.