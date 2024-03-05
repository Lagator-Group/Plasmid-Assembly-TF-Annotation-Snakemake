# %%
import pandas as pd
import os

# %%
plasmid_df = pd.read_csv('plasmid_summary.csv',sep=',')

# %%
n=0
try:
    os.mkdir('deepTFactor_hits')
except:
    pass

for locid_df in plasmid_df['DeepTFactor Locus_ID']:
    locid_list = locid_df.split(',')
    plasmid = plasmid_df.iloc[n]['Plasmid']
    print(plasmid)
    deepTFactor_hits = 'deepTFactor_hits/' + plasmid
    for locid in locid_list:
        locid = locid.replace('[','').replace(']','').replace("'","").replace(' ','')
        locid = '>'+locid
        with open('prokka/' + plasmid + '/'+plasmid+'.ffn','r') as f:
            data = f.read()
            start = data.find(locid)
            end = data.find('>', start+1)
            result = data[start:end]
            with open(deepTFactor_hits + '.fasta', 'a') as f:
                f.write(result)
            f.close()
    n+=1


