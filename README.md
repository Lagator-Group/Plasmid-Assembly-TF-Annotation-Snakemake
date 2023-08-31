# Sequence Assembly and Annotation

Uses Snakemake pipeline for sequence alignment and annotation. Needs Snakemake environments to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Complete Pipeline Instructions

First run will take the longest as it will need to install the various environments. All subsequent runs SHOULD be faster as the environments will already have been installed, and will just need to be run.

1. Create folder ```data``` and copy sequence files to it.
2. Rename the sequence files to follow the following naming convention:
```
sample_1.fastq #forward read
sample_2.fastq #reverse read
sample.fastq #long read
``` 
3. Open ```snakefile``` in a text editor and edit the value for sample. If you have many samples, it's best to run each pipeline in sequence instead of parrallel as it would seem to cause conflicts whilst being run. Further testing still required to understand how it all behaves.
4. In the ```snakemake``` conda environment, run ```snakemake --use-conda --cores all --conda-frontend conda```.
5. Results will be output in ```results``` folder.
6. If you spot an error in your script(s), it is best to let it fail naturally rather than cancelling it whilst it's running.

## Branch Menu
- Main: Hybrid assembly with unicycler, separation of plasmid contigs, swissprot and prokka annotation in ```results/annotation/```.
- Long assembly: Same as above, but long-read assembly only using flye.
- Short assembly Spades: Same as above, but short-read assembly only using Shovill and Spades.
- Short assembly Velvet: Same as above, but short-read assembly only using Shovill and Velvet.
- No Assembly: Same as above, but no assembly. Place assembled sequences in ```data/{sample}.fasta```.