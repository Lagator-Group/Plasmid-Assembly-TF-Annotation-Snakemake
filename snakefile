SAMPLES=['1284','7764','9246']

rule all:
    input:
        expand('results/prokka_plasmid/{sample}',sample=SAMPLES)

rule filtlong:
    input:
        'data/{sample}.fastq'
    output:
        temp('data/filtlong/{sample}.fastq.gz')
    conda:
        'env/filtlong.yml'
    shell:
        'filtlong --min_length 1000 --keep_percent 95 --target_bases 500000000 {input} | gzip > {output}'

#for unicycler rules, comment out the unnecessary rules depending on whether or not you want to run hybrid, long or short de novo assembly.

rule unicycler: #use to only assemble long reads
    input:
        'data/filtlong/{sample}.fastq.gz'
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    shell:
        'unicycler -l {input} -o {output}'
        
rule abricate:
    input:
        'results/unicycler/{sample}'
    output:
        'results/abricate_plasmid/{sample}.tab'
    conda:
        'env/abricate.yml'
    shell:
        'abricate -db plasmidfinder {input}/assembly.fasta > {output}'

rule contig_plasmid:
    input:
        'results/abricate_plasmid/{sample}.tab'
    output:
        'results/contigs_plasmid/{sample}.fasta'
    script:
        'bin/contig_plasmid.py'

rule prokka:
    input:
        'results/contigs_plasmid/{sample}.fasta'
    output:
        directory('results/prokka_plasmid/{sample}')
    conda:
        'env/prokka.yml'
    shell:
        'prokka --outdir {output} --prefix {wildcards.sample} {input}'