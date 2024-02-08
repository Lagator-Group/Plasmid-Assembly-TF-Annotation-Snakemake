conda activate flye

for fastq in fastq/*.fastq ; do
    sample=$(basename $fastq .fastq)
    mkdir flye
    mkdir flye/$sample
    flye --nano-raw $fastq -o flye/$sample -t 8
    mkdir fasta_wgs
    cp flye/$sample/assembly.fasta fasta_wgs/$sample\.fasta
done