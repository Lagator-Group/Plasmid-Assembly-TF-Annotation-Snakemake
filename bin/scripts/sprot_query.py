# sprot_query.py
# queries UniProt for SwissProt data
# output is a tab-delimited file with three columns:
#    Locus Tag, Entry, and Gene Name
# if no hits are found, output is a single-row DataFrame with a message
#
# input: tab-delimited BLAST results file
# output: tab-delimited file

import pandas as pd
import time

# create a DataFrame with the columns we need
query_df = pd.DataFrame(columns=['Locus Tag', 'Entry'])

plasmid = snakemake.wildcards.sample
# read in the BLAST results
blastx_df = pd.read_csv(snakemake.input[0], delimiter='\t')
#Ouput file
sprot_out = snakemake.output[0]

'''plasmid = 'CP025890.1'
blastx_df = pd.read_csv(f'blastx_best/{plasmid}.tsv', delimiter='\t')
sprot_out = f'sprot/{plasmid}.tsv'
'''
# check if there were no BLAST results
if blastx_df['Locus Tag'][0] == 'No BLAST results':
    # create a DataFrame with the same columns and a single row for "No BLAST results"
    query_df = pd.DataFrame({'Locus Tag': ['No BLAST results'],
                             'Entry': ['No BLAST results'],
                             'Protein names': ['No BLAST results'],
                             'Gene Names': ['No BLAST results']})

# if there were BLAST results
else:
    # create lists to store locus tags and entries
    locus_list = []
    entry_list = []

    # iterate over the locus tags and entries in the BLAST results
    n=0
    for locus in blastx_df['Locus Tag']:
        locus_list.append(locus)
        entry = str(blastx_df['Entry'][n])
        entry = entry.replace('sp|', '')
        end = entry.find('|')
        entry = entry[:end]
        entry_list.append(entry)
        print(f'{locus} mapped to {entry}')
        n += 1

    # create a DataFrame with the locus tags and entries
    query_df = pd.DataFrame({'Locus Tag': locus_list, 'Entry': entry_list})
    # create an empty DataFrame with the columns we want from UniProt
    df = pd.DataFrame(columns=['Entry', 'Protein names', 'Gene Names'])
    # iterate over the UniProt entries and try to get the data
    for i in query_df['Entry']:
        url = 'https://rest.uniprot.org/uniprotkb/' + i + '.tsv'
        df = pd.concat([df, pd.read_csv(url, sep='\t', usecols=['Entry', 'Protein names', 'Gene Names'])])
        print(f'{i} done')
        
    query_df = query_df.merge(df, on='Entry', how='left')

query_df = query_df.drop_duplicates(subset=['Locus Tag'])
# write the results to a new file
with open(sprot_out, 'w', newline='') as f:
    query_df.to_csv(f, index=False, sep='\t')

