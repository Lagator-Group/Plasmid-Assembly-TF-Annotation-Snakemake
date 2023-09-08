# Sequence Assembly and Annotation

Uses Snakemake pipeline for sequence alignment and annotation. Needs Snakemake environment to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Complete Pipeline Instructions
### Assembly: Hybrid, Short or Long
First run will take the longest as it will need to install the various environments. All subsequent runs SHOULD be faster as the environments will already have been installed, and will just need to be run.

1. Create folder ```data``` and copy sequence files to it.
2. Rename the sequence files to follow the following naming convention:
```
sample_1.fastq #forward read, only relevant for short read
sample_2.fastq #reverse read, only relevant for short read
sample.fastq #long read, only relevant for long read
``` 
3. Run ```scritps/get_samples.py``` to generate sample list.
4. Copy sample list to ```config.yml```
5. If necessary, adjust the values in ```config.yml```, particularly the threads available.
6. To execute:
```
snakemake -s {desired_snakefile_pipeline} --cores {core_available i.e 8} --use-conda --conda-frontend conda
```
7. Results will be generated and stored in ```results/```

### Annotation Only
Same as above, but place pre-assembled sequences in ```data/``` in ```.fasta``` format.

### Pipeline Options (-s)
- ```snakefile_hybrid```: Runs hybrid assembly using Unicycler.
- ```snakefile_long```: Runs long assembly using Flye. Input format can be changed in ```config.yml```
- ```snakefile_short```: Runs short assembly using Shovill + Skesa. If OS is Linux (i.e. not WSL), recommended to change assembler to ```spades``` in ```config.yml```. 
- ```snakefile_plasmid_from_genome```: Extracts and annotates plasmid sequences from pre-assembled WGS.
- ```snakefile_plasmid_annotation```: Annotates pre-assembled plasmid sequences.
- ```snakefile_plasmid_TF```: Annotates pre-assembled plasmid sequences and isolates transcription factor-related annotations to ```results/annotation/TF_{sample}.tsv```.