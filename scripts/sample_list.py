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

with open('sample_list.txt','w') as f:
    f.write('SAMPLES='+str(final_fastq))
