import os

#gets fastq files in data directory
fastq_list=[]
for fastq in os.listdir('data'):
    if fastq.endswith('.fastq'):
        fastq_list.append(fastq)

pre_final_fastq=[]
for fastq in fastq_list:
    if fastq.endswith('_1.fastq'):
       fastq=fastq.replace('_1.fastq','')
       pre_final_fastq.append(fastq)
    elif fastq.endswith('_2.fastq'):
       fastq=fastq.replace('_2.fastq','')
       pre_final_fastq.append(fastq)
    elif fastq.endswith('.fastq'):
       fastq=fastq.replace('.fastq','')
       pre_final_fastq.append(fastq)

final_fastq=[]
for fastq in pre_final_fastq:
    if fastq in final_fastq:
        continue
    else:
        final_fastq.append(fastq)

with open('denovo_plasmid.txt','w') as f:
    f.write('\
#!/bin/bash --login\n\
#$ -cwd\n\
\n\
#$ -l mem512 #gives 32GB memory\n\
#$ -pe smp.pe 8 #allocates 8 cores\n\
\n\
conda activate snakemake\n')
    
for fastq in final_fastq:
    with open('denovo_plasmid.txt','a') as f:
        f.write('snakemake --use-conda --cores all --conda-frontend conda results/prokka_plasmid/'+fastq+'\n')

#converts CRLF to LF
import re

def convert_crlf_to_lf(file_path):
    with open(file_path, 'r', newline='') as file:
        content = file.read()

    updated_content = re.sub('\r\n', '\n', content)

    with open(file_path, 'w', newline='') as file:
        file.write(updated_content)

# Example usage
convert_crlf_to_lf('denovo_plasmid.txt')
