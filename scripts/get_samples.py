import os

_fasta_list = []
for fasta in os.listdir('data'):
    _fasta_list.append(fasta)

fasta_list=[]
for fasta in _fasta_list:
    fasta=fasta.replace('.fasta', '')
    fasta_list.append(fasta)

with open('sample_list','w') as f:
    f.write(str(fasta_list))