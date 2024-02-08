conda activate flye

for fastq in fastq/*.fastq ; do
    sample=$(basename $fastq .fastq)
    flye --nano-raw $fastq -o flye/$sample
    mkdir fasta_wgs
    mv flye/$sample/assembly.fasta fasta_wgs/$sample\.fasta
done