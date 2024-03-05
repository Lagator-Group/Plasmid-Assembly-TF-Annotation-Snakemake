import pandas as pd
import os

prokka = snakemake.input[0]

plasmid_file = snakemake.input[1]
plasmid_name = plasmid_file.replace('.fasta','').replace('fasta_plasmid/','')

plasmid_df = pd.read_csv('plasmid_summary.csv',sep=',')

deepTFactor_hits = snakemake.output[0]

locid_list = plasmid_df.loc[plasmid_df['Plasmid']==plasmid_name]['DeepTFactor Locus_ID'].iloc[0]
locid_list = locid_list.split(',')

for locid in locid_list:
    locid = locid.replace('[','').replace(']','').replace("'","").replace(' ','')
    locid = '>'+locid
    with open(prokka,'r') as f:
        data = f.read()
        start = data.find(locid)
        end = data.find('>', start+1)
        result = data[start:end]
        with open(deepTFactor_hits, 'a') as f:
            f.write(result)
        f.close()
