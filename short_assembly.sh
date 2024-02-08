
conda activate shovill

for fastq in fastq/*_1.fastq ; do
    sample=$(basename $fastq _1.fastq)
    shovill --R1 fastq/$sample\_1.fastq --R2 fastq/$sample\_2.fastq --outdir shovill/$sample --force
    mkdir fasta_wgs
    mv shovill/$sample/contigs.fa fasta_wgs/$sample\.fasta
done