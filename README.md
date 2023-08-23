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
3. Open ```snakefile``` in a text editor and edit the value for sample.
4. Depending on whether you want to assemble short reads, long reads or do a hybrid assembly, comment out the ```rule unicycler``` rules you do not need.
5. In the ```snakemake``` conda environment, run ```snakemake --use-conda --cores all --conda-frontend conda```.
6. Results will be output in ```results``` folder.
7. If you want to process multiple samples, run the following command in the terminal: ```snakemake --use-conda --cores all --conda-frontend conda results/prokka_plasmid/{sample1} ; results/prokka_plasmid/{sample2} ; etc```. This is so that if 1 sample fails for a reason, it won't affect the processing of the other samples.
8. If you spot an error in your script(s), it is best to let it fail naturally rather than cancelling it whilst it's running.