# %%
import pandas as pd

# %%
locus_tags = []
tsv= snakemake.input[0]

df=pd.read_csv(tsv,sep='\t',usecols=['locus_tag','product'])
_n=df.shape[0] #get number of rows in df
n=0

while n<_n: #will only run for number of rows in pd.df
    locus=str(df['locus_tag'][n]) #looks at reach row 'n' of pd.df
    product='>'+str(df['product'][n]) #looks at corresponding contig code
    
    if product == '>hypothetical protein':
        locus_tags.append(locus)

    n=n+1

# %%
fnn= snakemake.input[1]
hypothetical_proteins = snakemake.output[0]
with open(fnn, 'r') as f:
    data=f.read()
    for locus in locus_tags:
        start=data.find(locus)
        end=data.find('>', start+1)
        result=data[start:end]

        with open(hypothetical_proteins,'a') as f:
            f.write('>'+result)
