# Cross-references all sprot hits with known DNA- and RNA-binding Transcription Regulators

import pandas as pd
import sprot_entries

# Read the .tsv file into a DataFrame
sprot_df = pd.read_csv(snakemake.input[0], delimiter='\t')

entries = sprot_entries.nucleotide_synthesis

# checks if entry is in above list of known DNA- or RNA-binding Transcription Regulators

# handles empty BLAST results
if sprot_df['Entry'][0] == 'No BLAST results':
    tf_df = pd.DataFrame({'Locus Tag':['No BLAST results'],'Entry':['No BLAST results'],'Protein names':['No BLAST results'],'Gene Names':['No BLAST results']})

else:
    # create lists to store locus tags, entries, protein names, and gene names
    locus_tag_list = []
    entry_list = []
    protein_list = []
    gene_list = []

    n = 0

    for entry in sprot_df['Entry']:
        # iterate over each entry in the BLAST results
        print(f"Searching for {entry}")

        # if the current entry is in the list of known TFs and the locus tag has not been added yet,
        # add the locus tag, entry, protein name, and gene name to their respective lists
        if entry in entries and sprot_df['Locus Tag'][n] not in locus_tag_list:
            print(f"Found {entry}! Adding to lists.")
            locus_tag_list.append(sprot_df['Locus Tag'][n])
            entry_list.append(sprot_df['Entry'][n])
            protein_list.append(sprot_df['Protein names'][n])
            gene_list.append(sprot_df['Gene Names'][n])

        # otherwise do nothing
        else:
            print(f"{entry} not found or locus tag already added.")
        n += 1

    # create a dictionary with the lists as values
    tf_data = {
        'Locus Tag': locus_tag_list,
        'Entry': entry_list,
        'Protein names': protein_list,
        'Gene Names': gene_list
    }

    # create a pandas DataFrame from the dictionary
    tf_df = pd.DataFrame(tf_data)

with open(snakemake.output[0],'w',newline='') as f:
    tf_df.to_csv(f,index=False,sep='\t')
