conda activate unicycler

for fastq in fastq/*_1.fastq ; do
    sample=$(basename $fastq _1.fastq)
    unicycler -1 fastq/$sample\_1.fastq -2 fastq/$sample\_2.fastq -l fastq/$sample\.fastq -o unicycler/$sample
    mkdir fasta_wgs
    mv unicycler/$sample/assembly.fasta fasta_wgs/$sample\.fasta
done