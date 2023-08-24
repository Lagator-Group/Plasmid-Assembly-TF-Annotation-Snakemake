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
5. In the ```snakemake``` conda environment, run ```snakemake --use-conda --cores all --conda-frontend conda results/{sample1}```.
6. Results will be output in ```results``` folder.
7. If you want to process multiple samples, run the following command in the terminal: ```snakemake --use-conda --cores all --conda-frontend conda results/{sample1} ; snakemake --use-conda --cores all --conda-frontend conda results/{sample2} ; etc```. This is so that if 1 sample fails for a reason, it won't affect the processing of the other samples.
8. If you spot an error in your script(s), it is best to let it fail naturally rather than cancelling it whilst it's running.

## Preparing batch files
I've prepared scripts ```make_batch.py``` to simply the preparation of the batch files if you have many samples you want to analyse. Run the following command in the terminal:
```
python3 make_batch.py
```
All the necessary batch files will get written (in the right format), and file ```qsub_commands.txt``` will contain the command to queue them. Just paste the text into the command line to queue all your jobs.