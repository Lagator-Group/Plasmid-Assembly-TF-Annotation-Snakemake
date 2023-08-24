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

#writes batch file
def write_batch_file(sample):
    with open('assembly_'+sample+'.txt','w') as f:
        f.write('\
#!/bin/bash --login\n\
#$ -cwd\n\
#$ -l mem512 #gives 32GB memory\n\
#$ -pe smp.pe 8 #allocates 8 cores\n\
\n\
conda activate snakemake\n\
snakemake --use-conda --conda-frontend --cores all conda results/'+sample)

for fastq in final_fastq:
    write_batch_file(fastq)

import re

def convert_crlf_to_lf(file_path):
    with open(file_path, 'r', newline='') as file:
        content = file.read()

    updated_content = re.sub('\r\n', '\n', content)

    with open(file_path, 'w', newline='') as file:
        file.write(updated_content)

for file in os.listdir():
    if file.endswith('.txt'):
        convert_crlf_to_lf(file)


#prepares qsub commands for pasting into terminal
batch=[]
for file in os.listdir():
    if file.endswith('.txt'):
        batch.append('qsub '+file)

with open('qsub_commands.txt','w') as f:
    f.write(';'.join(batch))