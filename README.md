# Sequence Assembly and Annotation
Uses Snakemake pipeline for sequence alignment and annotation. Needs Snakemake environment to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Complete Pipeline Instructions

### Environments
All snakefiles can be run without pre-installing the necessary conda environments. The versions that do not require environment pre-installation have suffix `noenv`. This means that the environments will be installed within the `.snakemake` folder. These can be deleted manually after completion. If you wish to install the environments yourself, do so then run the snakefiles without the `noenv` suffix.

### 1. Get samples list (Optional)
Only needs to be run if you have too many samples to manually write in config file.
Will only get samples from folders further in pipeline (i.e. `fasta_plasmid` > `fasta_wgs` > `fastq`).
```
python3 bin/scripts/get_samples.py
```
Samples will be output to `sample_list.txt`. Copy sample list into config file.
### 2. Assembly: Hybrid, Short or Long
If your sequences are assembled, skip to next section. If you need to install conda environments, they can all be found in `bin/env` directory. Install with `conda env create -f path/to/env.yml`.
1. Create `fastq` directory.
2. Copy **ONLY** relevant fastq to folder (e.g. only copy `_1.fastq` and `_2.fastq` if assembling short-read sequences).
3. Run the appropriate assembly assembly script (e.g. `snakemake -s hybrid_assembly -c8 --use-conda --conda-frontend conda`).
4. Asssembled `.fasta` will be output to `fasta_wgs` directory.

### 3. Transcription Factor annotation
Ensure `config.yml` contains the proper sample names in the proper format.
```
samples: ['sample_1','sample_2']
```
If extracting plasmid TFs from WGS `.fasta`, ensure sequences are in `fasta_wgs` folder. If working with plasmid sequences, place `.fasta` in `fasta_plasmid` then run:
```
snakemake -s snakefile_plasmid_tf -c8 --use-conda --conda-frontend-conda
```

### 4. DNA-binding motif prediction
**Bug**: *This step does not currently work on the UoM CSF, but should work on local machines relatively quickly.*
Ensure you have installed the `deeptfactor` conda env from `bin/env/deeptfactor.yml`.
Script requires `fasta_plasmid` directory to be populated already as it gets the plasmid names from there.
Run with `snakemake -s snakefile_deeptfactor -c8 --use-conda`.
Results will be output in `deeptfactor` directory.

A summary of the results will output to `plasmid_summary.csv`.