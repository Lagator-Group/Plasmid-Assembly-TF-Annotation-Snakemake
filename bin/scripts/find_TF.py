import pandas as pd

# Read the .tsv file into a DataFrame
sprot_df = pd.read_csv(snakemake.input[0], delimiter='\t')

tf_names=['transcription','repressor','activator','regulator']

locus_tag_list=[]
entry_list=[]
protein_list=[]
gene_list=[]
n=0
_n=sprot_df.shape[0]

while n<_n:
    for name in tf_names:
        if name in sprot_df.iloc[n]['Protein names'] and sprot_df.iloc[n]['Protein names'] not in protein_list:
            locus_tag_list.append(sprot_df.iloc[n]['Locus Tag'])
            entry_list.append(sprot_df.iloc[n]['Entry'])
            protein_list.append(sprot_df.iloc[n]['Protein names'])
            gene_list.append(sprot_df.iloc[n]['Gene Names'])
        else:
            pass
    n+=1
tf_data={'Locus Tag':locus_tag_list,'Entry':entry_list,'Protein names':protein_list,'Gene Names':gene_list}
tf_df=pd.DataFrame(tf_data)

with open(snakemake.output[0],'w') as f:
    tf_df.to_csv(f,index=False,sep='\t')