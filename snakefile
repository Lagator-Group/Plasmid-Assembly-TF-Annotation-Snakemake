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
    
rule get_hypothetical:
    input:
        tsv='results/prokka_plasmid/{sample}/{sample}.tsv',
        fnn='results/prokka_plasmid/{sample}/{sample}.fnn'
    output:
        temp('results/annotation/{sample}_hypothetical.fasta')
    script:
        'scripts/get_hypothetical.py'

rule uniprot_sprot:
    input:

    output:
        'bin/swissprot/uniprot_sprot.fasta'
    params:
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
    shell:
        'wget {params} -P bin/swissprot/ && '
        'gunzip {output}.gz'
'''
rule: makeblastdb:
    input:
        uniprot_fasta='bin/swissprot/uniprot_sprot.fasta'
    output:
    conda:
    log:
    shell:
    'makeblastdb -dbtype prot -in {output.uniprot_fasta} &&'
    'blastx -query results/annotation/7764_hypothetical.fasta -db bin/swissprot'
    ' -out results/blast/7764.tsv -outfmt 6 -evalue 10-max_hsps 1 -num_threads 4'
'''