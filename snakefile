SAMPLES=['1284']

rule all:
    input:
        expand('data/{sample}.fastq.gz', sample=SAMPLES)

rule filtlong:
    input:
        'data/{sample}.fastq'
    output:
        'data/{sample}.fastq.gz'
    conda:
        'env/filtlong.yml'
    log:
        'log/filtlong/{sample}.log'
    shell:
        '(filtlong --min_length 1000 --keep_percent 95 --target_bases 500000000 {input} | gzip > {output}) > {log} 2>&1'

#for unicycler rules, comment out the unnecessary rules depending on whether or not you want to run hybrid, long or short de novo assembly.
rule unicycler: #use for hybrid assembly
    input:
        short_F='data/{sample}_1.fastq',
        short_R='data/{sample}_2.fastq',
        _long='data/{sample}.fastq.gz'
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    shell:
        'unicycler -1 {input.short_F} -2 {input.short_R} -l {input._long} -o {output}'


'''
rule unicycler: #use to only assembly short reads
    input:
        short_F='data/{sample}_1.fastq',
        short_R='data/{sample}_2.fastq',
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    shell:
        'unicycler -1 {input.short_F} -2 {input.short_R} -o {output}'
        
rule unicycler: #use to only assemble long reads
    input:
        _long='data/{sample}.fastq.gz'
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    shell:
        'unicycler -l {input._long} -o {output}'
'''
        
rule abricate:
    input:
        'results/unicycler/{sample}/assembly.fasta'
    output:
        'results/abricate_plasmid/{sample}.tab'
    conda:
        'env/abricate.yml'
    log:
        'log/abricate_plasmid/{sample}.log'
    shell:
        '(abricate -db plasmidfinder {input} > {output}) > {log} 2>&1'

rule contig_plasmid:
    input:
        'results/abricate_plasmid/{sample}.tab'
    output:
        'results/contigs_plasmid/{sample}.fasta'
    log:
        'log/contig_plasmid/{sample}.log'
    script:
        'bin/contig_plasmid.py'

rule prokka:
    input:
        'results/contigs_plasmid/{sample}.fasta'
    output:
        directory('results/prokka_plasmid/{sample}')
    conda:
        'env/prokka.yml'
    log:
        'log/prokka_plasmid/{sample}.log'
    shell:
        '(prokka --outdir {output} --prefix {sample} {input}) > {log} 2>&1'