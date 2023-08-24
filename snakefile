SAMPLES=['1284','7764','9246']

rule all:
    input:
        expand('results/{sample}',sample=SAMPLES)

#for unicycler rules, comment out the unnecessary rules depending on whether or not you want to run hybrid, long or short de novo assembly.
'''
rule unicycler: #hybrid assembly
    input:
        long1='data/{sample}.fastq',
        short1='data/{sample}_1.fastq',
        short2='data/{sample}_2.fastq'
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    shell:
        '(unicycler -1 {input.short1} -2 {input.short2} -l {input.long1} -o {output}) > {log}'

rule unicycler: #short assembly
    input:
        short1='data/{sample}_1.fastq',
        short2='data/{sample}_2.fastq'
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    shell:
        '(unicycler -1 {input.short1} -2 {input.short2} -o {output}) > {log}'
'''
rule unicycler: #long assembly
    input:
        'data/{sample}.fastq'
    output:
        directory('results/unicycler/{sample}')
    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    shell:
        '(unicycler -l {input} -o {output}) > {log}'

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
    log:
        'log/contig_plasmid/{sample}.log'
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

rule end:
    input:
        'results/prokka_plasmid/{sample}'
    output:
        'results/{sample}'
    shell:
        'echo "done" > {output}'