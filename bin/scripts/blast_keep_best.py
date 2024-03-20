
import pandas as pd

# Read the .tsv file into a DataFrame
blastx_df = pd.read_csv(snakemake.input[0],header=None,delimiter='\t')

locus_list = blastx_df[0].to_list()
entry_list = blastx_df[1].to_list()

best_locus = []
best_entry = []

n=0
for locus in locus_list:
    if locus in best_locus:
        pass
    else:
        best_locus.append(locus)
        best_entry.append(entry_list[n])
        print(locus+' '+entry_list[n])
    n+=1

best_df = pd.DataFrame({'Locus Tag':best_locus,'Entry':best_entry})

with open(snakemake.output[0],'w',newline='') as f:
    best_df.to_csv(f,index=False,sep='\t')