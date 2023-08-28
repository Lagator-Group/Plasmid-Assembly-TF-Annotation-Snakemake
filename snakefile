configfile: 'config.yml'

rule all:
    input:
        'results/swissprot/1284_names.tsv'

rule unicycler: 
    input:
        long_='data/{sample}.fastq'
    output:
        folder=directory('results/unicycler/{sample}'),
        fasta='results/unicycler/{sample}/assembly.fasta'

    conda:
        'env/unicycler.yml'
    log:
        'log/unicycler/{sample}.log'
    params:
        keep='0',
        mode='conservative'
    threads: 
        config['threads']
    shell:
        'unicycler -l {input.long_} -o {output.folder} --keep {params.keep} -t {threads} --mode {params.mode}'
'''
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
'''
rule abricate:
    input:
        'results/unicycler/{sample}/assembly.fasta'
    output:
        'results/abricate_plasmid/{sample}.tab'
    conda:
        'env/abricate.yml'
    log:
        'log/abricate/{sample}.log'
    shell:
        'abricate -db plasmidfinder {input} > {output}'

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
        folder=directory('results/prokka_plasmid/{sample}'),
        tsv='results/prokka_plasmid/{sample}/{sample}.tsv',
        ffn='results/prokka_plasmid/{sample}/{sample}.ffn'

    conda:
        'env/prokka.yml'
    log:
        'log/prokka/{sample}.log'
    shell:
        'prokka --outdir {output.folder} --prefix {wildcards.sample} {input} --force'
    
rule get_hypothetical:
    input:
        tsv='results/prokka_plasmid/{sample}/{sample}.tsv',
        ffn='results/prokka_plasmid/{sample}/{sample}.ffn'
    output:
        'results/annotation/{sample}_hypothetical.fasta'
    script:
        'scripts/get_hypothetical.py'

rule wget_uniprot:
    input:

    output:
        uniprot_fasta='bin/swissprot/uniprot_sprot.fasta',
    params:
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
    log:
        'log/wget_uniprot.log'
    shell:
        'wget {params} -P bin/swissprot && '
        'gunzip {output.uniprot_fasta}.gz'

rule makeblastdb:
    input:
        uniprot_fasta='bin/swissprot/uniprot_sprot.fasta'
    output:
        phr='bin/swissprot/swissprot.phr',
        pin='bin/swissprot/swissprot.pin',
        psq='bin/swissprot/swissprot.psq'
    conda:
        'env/blast.yml'
    log:
        'log/makeblastdb.log'
    params:
        db='bin/swissprot/swissprot'
    shell:
        'makeblastdb -dbtype prot -in {input.uniprot_fasta} -out {params.db}'

rule blastx:
    input:
        hypothetical='results/annotation/{sample}_hypothetical.fasta',
        phr='bin/swissprot/swissprot.phr',
        pin='bin/swissprot/swissprot.pin',
        psq='bin/swissprot/swissprot.psq'
    output:
        'results/swissprot/{sample}.tsv'   
    conda:
        'env/blast.yml'
    log:
        'log/blastx/{sample}.log'
    params:
        db='bin/swissprot/swissprot',
        outfmt='6',
        evalue='0.1',
        max_hsps='1'
    threads:
        config['threads']
    shell:
        'blastx -query {input.hypothetical} -db {params.db}'
        ' -out {output} -outfmt {params.outfmt} -evalue {params.evalue} -max_hsps {params.max_hsps} -num_threads {threads}'

rule uniprot_query:
    input:
        'results/swissprot/{sample}.tsv'
    output:
        'results/swissprot/{sample}_names.tsv'
    script:
        'scripts/uniprot_query.py'