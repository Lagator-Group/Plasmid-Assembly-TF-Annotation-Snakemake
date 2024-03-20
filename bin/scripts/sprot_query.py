import pandas as pd
import time

query_df=pd.DataFrame(columns=['Locus Tag','Entry'])
query_df

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

error = 0
error_max = 10
while error <= error_max:
    try:
        for i in query_df['Entry']:
            print(i)
            url = 'https://rest.uniprot.org/uniprotkb/'+i+'.tsv'
            df = pd.concat([df,pd.read_csv(url, sep='\t',usecols=['Entry','Protein names','Gene Names'])])
    except:
        time.sleep(60)
        error += 1

query_df=query_df.merge(df,on='Entry',how='left')

with open(snakemake.output[0],'w') as f:
    query_df.to_csv(f,index=False,sep='\t')


