# Sequence Assembly and Annotation
Uses Snakemake pipeline for sequence alignment and annotation. Needs Snakemake environment to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Description
These scripts are designed to assemble NGS/ONT raw with [Shovill](https://github.com/tseemann/shovill), [Flye](https://github.com/fenderglass/Flye), or [Unicyler](https://github.com/rrwick/Unicycler) for short-read, long-read or hybrid assembly, respectively. Pipeline will then separate plasmid contigs from WGS `.fasta`, then annotate the plasmid sequence. Plasmid contigs and AMR markers are identified using [Abricate](https://github.com/tseemann/abricate). Annotation is done by identifying all CDS with [Prokka](https://github.com/tseemann/prokka) and using BLAST against the Swissprot database. Highest ranking hits are kept and cross-references with all known DNA- and RNA-binding Transcription Regulators in the Swissprot database. All CDS are also process with [DeepTFactor](https://bitbucket.org/kaistsystemsbiology/deeptfactor/src/master/). All data is then summarised in `plasmid_summary.csv`.

## Complete Pipeline Instructions

### Environments
All snakefiles can be run without pre-installing the necessary conda environments. The versions that do not require environment pre-installation have suffix `noenv`. This means that the environments will be installed within the `.snakemake` folder, which can be deleted manually after completion. If you wish to install the environments yourself, do so then run the snakefiles without the `noenv` suffix.

### 1. Config (Optional)
Adjust the thread number to what your machine is capable of.

Run `python3 bin/scripts/get_samples.py` to automatically overwrite `config.yml` with proper sample names.
Will only get samples from folders further in pipeline (i.e. `fasta_plasmid` > `fasta_wgs` > `fastq`).

#### Folder Description
`fastq`: Put raw NGS/ONT reads here.
`fasta_wgs`: Put assembled WGS here.
`fasta_plasmid`: Put assembled PLASMID sequences here.

### 2. Assembly: Hybrid, Short or Long
If your sequences are assembled, skip to next section.
1. Copy **ONLY** relevant fastq to folder (e.g. only copy `_1.fastq` and `_2.fastq` if assembling short-read sequences).
2. Run the appropriate assembly assembly script (e.g. `snakemake -s hybrid_assembly -c8 --use-conda --conda-frontend conda`).
3. Asssembled `.fasta` will be output to `fasta_wgs` directory.

### 3. Transcription Factor annotation
If extracting plasmid TFs from WGS `.fasta`, ensure sequences are in `fasta_wgs` folder. If working with plasmid sequences, place `.fasta` in `fasta_plasmid` then run:
```
snakemake -s snakefile_plasmid_tf -c8 --use-conda --conda-frontend conda
```

### 4. DNA-binding motif prediction
**Bug**: *This step does not currently work on the UoM CSF, but should work on local machines relatively quickly. Just make sure to download all the folders that have been generated as they are needed for this step.*
Ensure you have installed the `deeptfactor` conda env from `bin/env/deeptfactor.yml`.
Run with `snakemake -s snakefile_deeptfactor -c8 --use-conda --conda-frontend conda`.
Results will be output in `deeptfactor` directory.

A summary of the results will output to `plasmid_summary.csv`.