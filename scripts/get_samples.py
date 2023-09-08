import os

file_list = []
for file in os.listdir('data'):
    file_list.append(file)

_fasta_list=[]
_fastq_list=[]

for file in file_list:
    if file.endswith('.fasta'):
        _fasta_list.append(file)
    elif file.endswith('.fastq'):
        _fastq_list.append(file)

if len(_fasta_list) > 0 and len(_fastq_list) > 0:
    print('Found both .fasta and .fastq files, please keep 1 data format in the folder')
    quit()

fasta_list=[]
for fasta in _fasta_list:
    fasta=fasta.replace('.fasta', '')
    fasta_list.append(fasta)

long_list=[]
_1_list=[]
_2_list=[]
for fastq in _fastq_list:
    if fastq.endswith('_1.fastq'):
        _1_list.append(fastq)
    elif fastq.endswith('_2.fastq'):
        _2_list.append(fastq)
    else:
        long_list.append(fastq)

fastq_list=[]
if  (len(_1_list) == len(_2_list) and len(_1_list) == 2*len(long_list)) or \
    (len(long_list) > 0 and len(_1_list) == 0 and len(_2_list) == 0):
    for fastq in long_list:
        fastq.replace('.fastq', '')
        fastq_list.append(fastq)
elif len(_1_list) == len(_2_list) and len(long_list) == 0:
    for fastq in _1_list:
        fastq.replace('_1.fastq', '')
        fastq_list.append(fastq)
else:
    print('Something went wrong, please check your files')
    quit()

if len(fastq_list) > 0 and len(fasta_list) == 0:
    sample_list = fastq_list
elif len(fasta_list) > 0 and len(fastq_list) == 0:
    sample_list = fasta_list
else:
    print('Something went wrong, please check your files')
    quit()

with open('sample_list','w') as f:
    f.write(str(sample_list))