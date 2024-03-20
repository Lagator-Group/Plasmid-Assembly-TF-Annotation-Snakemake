import pandas as pd
import time

query_df=pd.DataFrame(columns=['Locus Tag','Entry'])
query_df

plasmid = snakemake.input[0]

with open(snakemake.input[0]) as f:
    for line in f:
        locus_tag=line[:14]
        start=line.find('sp|')
        end=line.find('|',start+3)
        ID=line[start:end]
        ID=ID.replace('sp|','')

        new_row={'Locus Tag':locus_tag,'Entry':ID}
        query_df=pd.concat([query_df,pd.DataFrame([new_row])])

query_df.reset_index(drop=True,inplace=True)

df=pd.DataFrame(columns=['Entry','Protein names','Gene Names'])

error_max = 10
for i in query_df['Entry']:
    error = 0
    if error == error_max:
        print(f'Not able to process entry {i} on plasmid {plasmid} after {error} tries')
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


