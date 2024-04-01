import pandas as pd
import os

# read in Prokka file
prokka = snakemake.input[0]
print(f'prokka file: {prokka}')

# read in plasmid file name
plasmid_file = snakemake.input[1]
plasmid_name = plasmid_file.replace('.fasta','').replace('fasta_plasmid/','')
print(f'plasmid name: {plasmid_name}')

# read in plasmid summary file
plasmid_df = pd.read_csv('plasmid_summary.csv',sep=',')
print(f'plasmid_df shape: {plasmid_df.shape}')

# read in DeepTFactor hits
deepTFactor_hits = snakemake.output[0]
print(f'deepTFactor_hits: {deepTFactor_hits}')

# find locus IDs for the current plasmid
locid_list = plasmid_df.loc[plasmid_df['Plasmid']==plasmid_name]['DeepTFactor Locus_ID'].iloc[0]
print(f'locid_list: {locid_list}')
locid_list = locid_list.split(',')
print(f'locid_list len: {len(locid_list)}')

# iterate over locus IDs and extract corresponding sequences from Prokka file
for locid in locid_list:
    print(f'locid: {locid}')
    # remove brackets, quotes, and whitespace
    locid = locid.replace('[','').replace(']','').replace("'","").replace(' ','')
    print(f'locid (cleaned): {locid}')
    # add '>' to start of locus ID
    locid = '>'+locid
    print(f'locid (with >): {locid}')
    with open(prokka,'r') as f:
        # read in Prokka file
        data = f.read()
        print(f'prokka file size: {len(data)}')
        # find start and end indices of locus ID in Prokka file
        start = data.find(locid)
        end = data.find('>', start+1)
        print(f'start: {start}, end: {end}')
        # extract sequence
        result = data[start:end]
        with open(deepTFactor_hits, 'a') as f:
            # write to output file
            f.write(result)
            print(f'wrote {len(result)} bytes to {deepTFactor_hits}')
        f.close()

