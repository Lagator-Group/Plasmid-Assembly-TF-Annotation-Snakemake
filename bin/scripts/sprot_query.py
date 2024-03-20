import pandas as pd
import time

query_df=pd.DataFrame(columns=['Locus Tag','Entry'])
query_df

plasmid = snakemake.input[0]

blastx_df = pd.read_csv(snakemake.input[0],delimiter='\t')

locus_list = []
entry_list = []
n=0
for locus in blastx_df['Locus Tag']:
    locus_list.append(locus)
    entry = blastx_df['Entry'][n]
    entry = entry.replace('sp|','')
    end = entry.find('|')
    entry = entry[:end]
    entry_list.append(entry)
    n+=1

query_df = pd.DataFrame({'Locus Tag':locus_list,'Entry':entry_list})
df=pd.DataFrame(columns=['Entry','Protein names','Gene Names'])

error_max = 10
for i in query_df['Entry']:
    error = 0
    if error == error_max:
        print(f'Not able to process entry {i} on {plasmid} after {error} tries')
        quit()
    else:
        try:
            url = 'https://rest.uniprot.org/uniprotkb/'+i+'.tsv'
            df = pd.concat([df,pd.read_csv(url, sep='\t',usecols=['Entry','Protein names','Gene Names'])])
            print(f'{i} done')
        except:
            print(f'error with query {i}, trying again in 60 seconds')
            time.sleep(60)
            error += 1

query_df=query_df.merge(df,on='Entry',how='left')

with open(snakemake.output[0],'w') as f:
    query_df.to_csv(f,index=False,sep='\t')


