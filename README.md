# Sequence Alignment and Annotation
## Environment Preparation
All scripts can only be run within Linux. All scripts were developed in WSL.
Instructions on how to install WSL can be found [here](https://learn.microsoft.com/en-us/windows/wsl/install).

It is strongly recommended to use conda for each package. Several bash scripts in this repo assume you have conda installed. A setup.sh script is included to create all the necessary environments and install the appropriate packages.
Instructions to install Anaconda can be found [here](https://gist.github.com/kauffmanes/5e74916617f9993bc3479f401dfec7da).

To create all the conda envs and install all the relevant packages run ```bash -i seq_tools/ngs_setup.sh``` or ```bash -i seq_tools/rnaseq_setup.sh``` and press "y" when prompted.

Before starting, open the ```config.ini``` file and change the values to match with your system settings

To remove the environments, run either ```bash -i ngs_uninstall.sh``` or ```bash -i rnaseq_uninstall.sh```

## Complete Pipeline Instructions
### Assemble and Isolate Plasmid Contigs
Copy the folder ```seq_tools``` to the folder containing the sequence files and run the following:
```
bash -i seq_tools/assembly_plasmid.sh
```

### RNAseq
Download reference genome from NCBI and copy `GCF###.fna` and `genome.gtf` files to directory containing `.fastq` files.
Copy the folder ```seq_tools``` to the folder containing the sequence files and run the following:
```
bash -i seq_tools/rnaseq.sh
```