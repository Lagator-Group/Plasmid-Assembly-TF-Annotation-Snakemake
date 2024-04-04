# Sequence Assembly and Annotation
Uses [Snakemake](https://github.com/snakemake/snakemake) pipeline for sequence alignment and annotation. Needs Snakemake environment to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Description
These scripts are designed to assemble NGS/ONT raw with [Shovill](https://github.com/tseemann/shovill), [Flye](https://github.com/fenderglass/Flye), or [Unicyler](https://github.com/rrwick/Unicycler) for short-read, long-read or hybrid assembly, respectively. Pipeline will then separate plasmid contigs from WGS `.fasta`, then annotate the plasmid sequence. Plasmid contigs and AMR markers are identified using [Abricate](https://github.com/tseemann/abricate). Annotation is done by identifying all CDS with [Prokka](https://github.com/tseemann/prokka) and using BLAST against the Swissprot database. Highest ranking hits are kept and cross-references with all known DNA- and RNA-binding Transcription Regulators in the Swissprot database. All CDS are also processed with [DeepTFactor](https://bitbucket.org/kaistsystemsbiology/deeptfactor/src/master/). All data is then summarised in `plasmid_summary.csv`.

## Complete Pipeline Instructions
All snakemake pipeline files (A.K.A. snakefiles) are prefixed with `sfile_`.

### Environments
All snakefiles can be run without pre-installing the necessary conda environments. The versions that do not require environment pre-installation have suffix `_noenv`. This means that the environments will be installed within the `.snakemake` folder, which can be deleted manually after completion. If you wish to install the environments yourself, do so then run the snakefiles without the `_noenv` suffix.

### Swissprot Database
`sfile_plasmid_tf` will download and set-up the Swissprot database in `bin/swissprot` the first time it is launched. It will not automatically do so on subsequent runs. If you wish to update the database, you can delete `bin/swissprot` before running the pipeline and it will automatically download and install the latest version.

### 1. Config
Adjust `config.yml` thread number to what your machine is capable of.

Run to automatically overwrite `config.yml` with proper sample names.
```
python3 bin/scripts/get_samples.py
```
Will only get samples from folders further in pipeline (i.e. `fasta_plasmid` > `fasta_wgs` > `fastq`).

#### Folder Description
`fastq`: Put raw NGS/ONT reads here.

`fasta_wgs`: Put assembled WGS here.

`fasta_plasmid`: Put assembled PLASMID sequences here.

### 2. Assembly: Hybrid, Short or Long
If your sequences are assembled, skip to next section.
1. Copy **ONLY** relevant fastq to folder (e.g. only copy `_1.fastq` and `_2.fastq` if assembling short-read sequences).
2. Run the appropriate assembly assembly script 
```
snakemake -s sfile_hybrid_assembly -c8 --use-conda --conda-frontend conda
snakemake -s sfile_short_assembly -c8 --use-conda --conda-frontend conda
snakemake -s sfile_long_assembly -c8 --use-conda --conda-frontend conda
```
3. Asssembled `.fasta` will be output to `fasta_wgs` directory.

### 3. Transcription Factor annotation
If extracting plasmid TFs from WGS `.fasta`, ensure sequences are in `fasta_wgs` folder. If working with plasmid sequences, place `.fasta` in `fasta_plasmid` then run:
```
snakemake -s sfile_plasmid_tf -c8 --use-conda --conda-frontend conda
```

#### Alternate identification
If you wish to identify different category of genes, make a list of all known entries from Uniprot and replace the entry names in `bin/scripts/sprot_entries.py` accordingly.

### 4. DNA-binding motif prediction
**Bug**: *This step does not currently work on the UoM CSF, but should work on local machines relatively quickly. Just make sure to download all the folders that have been generated as they are needed for this step.*

Ensure you have installed the `deeptfactor` conda env from `bin/env/deeptfactor.yml` if not using the `_noenv` sfile.
```
snakemake -s sfile_deeptfactor -c8 --use-conda --conda-frontend conda
```

A summary of the results will output to `plasmid_summary.csv`.

**NOTE**: deepTFactor cannot handle stop codons (*) in the `prokka.faa` CDS. Script will automatically make a copy without those. When surveying results, just quickly check that your finding does not have a stop codon by looking at the original `prokka` folder.

### Known Issues
**The environment for deepTFactor so run will not install properly on the UoM CSF.**

Download all the updated/new directories and run `sfile_deeptfactor` locally.


**`rule sprot_query` in `sfile_plasmid_tf` will crash on the CSF. This appears to be caused by a UoM networking issue.**

Adjust `rule all` in `sfile_plasmid_tf` to the following:
```
rule all:
    input:
        expand('abricate_amr/{sample}.tab',sample=config['samples']),
        expand('blastx/{sample}.tsv',sample=config['samples'])
```
This will run the offline (and resource-intensive) sections of the pipeline on the CSF. Download the `abricate`, `prokka` and `blastx` directories and run the pipeline again offline. It *should* be relatively quick. If you need to pause the pipeline at any stage, press `CTRL-Z` to halt the pipeline. Then run the following when ready to resume:
```
snakemake -s sfile_plasmid_tf --unlock
snakemake -s sfile_plasmid_tf -c8 --use-conda --rerun-incomplete
```