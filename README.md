### Current plans for improvement
1. To include DNA-binding protein prediction using deeptfactor. Necessary scripts, environments and instructions currently found in `bin/scripts/deeptfactor`.

# Sequence Assembly and Annotation
Uses Snakemake pipeline for sequence alignment and annotation. Needs Snakemake environment to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Complete Pipeline Instructions
### Assembly: Hybrid, Short or Long
If your sequences are assembled, skip to next section. If you need to install conda environments, they can all be found in `bin/env` directory. Install with `conda env create -f path/to/env.yml`.
1. Create `fastq` directory.
2. Copy **ONLY** relevant fastq to folder (e.g. only copy `_1.fastq` and `_2.fastq` if assembling short-read sequences).
3. Run the appropriate assembly `.sh` with `bash -i {*}_assembly.sh`.
4. Asssembled `.fasta` will be output to `fasta_wgs` directory.

### Transcription Factor annotation
If extracting plasmid TFs from WGS `.fasta`, run:
```
snakemake -s snakefile_tf_from_wgs -c8 --use-conda 
```
If plasmid `.fasta` already isolated, place sequences in ```fasta_plasmid``` directory and run:
```
snakemake -s snakefile_tf_from_plasmid -c8 --use-conda
```