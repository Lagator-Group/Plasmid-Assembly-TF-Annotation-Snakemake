SAMPLES=['1284','7764','9246']

rule all:
    input:
        expand('results/prokka_plasmid/{sample}',sample=SAMPLES)

rule unicycler: 
    input:
        long_='data/{sample}.fastq',
        short_1='data/{sample}_1.fastq',
        short_2='data/{sample}_2.fastq'
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    shell:
        '(unicycler -1 {input.short_1} -2 {input.short_2} -l {input.long_} -o {output}) > {log}'

rule abricate:
    input:
        'results/unicycler/{sample}'
    output:
        'results/abricate_plasmid/{sample}.tab'
    conda:
        'env/abricate.yml'
    log:
        'log/abricate/{sample}.log'
    shell:
        '(abricate -db plasmidfinder {input}/assembly.fasta > {output}) > {log}'

rule contig_plasmid:
    input:
        'results/abricate_plasmid/{sample}.tab'
    output:
        'results/contigs_plasmid/{sample}.fasta'
    script:
        'scripts/contig_plasmid.py'

rule prokka:
    input:
        'results/contigs_plasmid/{sample}.fasta'
    output:
        directory('results/prokka_plasmid/{sample}')
    conda:
        'env/prokka.yml'
    log:
        'log/prokka/{sample}.log'
    shell:
        '(prokka --outdir {output} --prefix {wildcards.sample} {input}) > {log}'