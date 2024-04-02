# tabulate_all.py
# Summarise all plasmid assemblies and annotation results in a single CSV
# Assumes that all files are in the right format
# Assemblies and annotations are in subdirectories:
#   fasta_plasmid/
#   fasta_wgs/
#   blastx/
#   abricate/
#   prokka/
#   deepTFactor/
# Output is plasmid_summary.csv


import pandas as pd
import os

# create a list of plasmid names from the fasta files
plasmid_list = []
for plasmid in os.listdir('fasta_plasmid'):
    plasmid = plasmid.replace('.fasta','')
    plasmid_list.append(plasmid)
    print('added plasmid: ' + plasmid)

# create a dataframe with plasmid names as rows and a single column
df=pd.DataFrame(columns=['Plasmid'])
df['Plasmid'] = plasmid_list
print('created dataframe with plasmid names')

# add column with plasmid length
plasmid_length = []
# iterate over plasmid names
for plasmid in df['Plasmid']:
    # open fasta file with plasmid name and read lines
    with open('fasta_plasmid/'+plasmid+'.fasta') as f:
        data = f.readlines()
        print('read fasta file for plasmid: ' + plasmid)
    # get length of second line (sequence)
    length = len(data[1])
    # add length to list
    plasmid_length.append(length)
    print('added length: ' + str(length) + ' to list')

# add plasmid length column to dataframe
df['Plasmid Length'] = plasmid_length
print('added plasmid length column to dataframe')

# lists to store all resistance and gene for each plasmid
resistance_all = []
gene_all = []
# list to store number of resistances for each plasmid
resistance_count = []

# iterate over plasmid names
for plasmid in df['Plasmid']:
    # print plasmid name for feedback
    print('current plasmid: ' + plasmid)
    # lists to store resistance and gene for current plasmid
    resistance_list = []
    gene_list = []
    
    # read amr_df for current plasmid
    amr_df = pd.read_csv('abricate_amr/'+plasmid+'.tab',sep='\t')
    
    # if no amr_df for current plasmid, add None to lists
    if len(amr_df) == 0:
        gene_list.append('None')
        resistance_list.append('None')
        resistance_count.append(0)
        print('no amr_df for plasmid: ' + plasmid)
    
    # if there is an amr_df for current plasmid
    else:
        # add all gene names to current lists
        for gene in amr_df['GENE']:
            gene_list.append(gene)
            print('added gene: ' + gene + ' to list')
        # add all resistance names to current lists
        for resistance in amr_df['RESISTANCE']:
            resistance_list.append(resistance)
            print('added resistance: ' + resistance + ' to list')
        # add number of resistances to current list
        resistance_count.append(len(resistance_list))
        print('added resistance count: ' + str(len(resistance_list)) + ' to list')

    # add current resistance and gene lists to all lists
    resistance_all.append(resistance_list)
    gene_all.append(gene_list)

# add resistance and gene lists as columns to dataframe
df['Resistance'] = resistance_all
df['Resistance Gene'] = gene_all
df['Resistance Count'] = resistance_count
print('added resistance and gene lists as columns to dataframe')

# read sprot_TF output for each plasmid
tf_prot_all = []
tf_gene_all = []
tf_count = []
for plasmid in df['Plasmid']:
    print('current plasmid: ' + plasmid)
    # lists to store TF protein and gene for current plasmid
    tf_prot_list = []
    tf_gene_list = []
    # read sprot_TF output for current plasmid
    tf_df = pd.read_csv('sprot_TF/'+plasmid+'.tsv',sep='\t')
    # if no sprot_TF output for current plasmid, add None to lists
    if len(tf_df) == 0:
        print('no sprot_TF output for plasmid: ' + plasmid)
        tf_gene_list.append('None')
        tf_prot_list.append('None')
        tf_count.append(0)
    # if there is sprot_TF output for current plasmid
    else:
        print('found sprot_TF output for plasmid: ' + plasmid)
        # add all gene names to current lists
        for gene in tf_df['Gene Names']:
            tf_gene_list.append(gene)
            print('added gene: ' + gene + ' to list')
        # add all protein names to current lists
        for prot in tf_df['Protein names']:
            tf_prot_list.append(prot)
            print('added protein: ' + prot + ' to list')
        # add number of TFs to current list
        tf_count.append(len(tf_prot_list))
        print('added TF count: ' + str(len(tf_prot_list)) + ' to list')
    # add current lists to all lists
    tf_prot_all.append(tf_prot_list)
    tf_gene_all.append(tf_gene_list)

# add TF lists to dataframe and simplify gene names
print("Simplifying gene names...")
tf_gene_simplified = []
for gene_list in tf_gene_all:
    print("Processing gene list " + str(len(tf_gene_simplified)+1) + "...")
    simplified_list = []
    for gene in gene_list:
        gene=str(gene)
        simplified_list.append(gene.split(' ')[0])
    tf_gene_simplified.append(simplified_list)
    print("Added gene list " + str(len(tf_gene_simplified)) + ".")

df['TF Protein'] = tf_prot_all
df['TF Gene'] = tf_gene_all
df['TF Gene Simplified'] = tf_gene_simplified
df['TF Count'] = tf_count

# read in DeepTF output for each plasmid
print("Reading in DeepTF output...")
deeptf_all = []
deeptf_count = []
for plasmid in df['Plasmid']:
    print("Currently processing plasmid " + plasmid + "...")
    deeptf_list = []  # list to store locus IDs for current plasmid
    deeptf_df = pd.read_csv('deepTFactor/'+plasmid+'/prediction_result.txt',sep='\t')  # read in DeepTF output
    n=0
    for prediction in deeptf_df['prediction']:  # iterate through DeepTF predictions
        prediction = str(prediction)  # convert to string
        if prediction == 'True':  # if prediction is TRUE, extract locus ID
            locid = deeptf_df['sequence_ID'][n]
            deeptf_list.append(locid)
        else:
            pass
        n+=1
    deeptf_count.append(len(deeptf_list))  # add number of locus IDs to current list
    deeptf_all.append(deeptf_list)  # add current list to all lists
    print("Added DeepTF list for plasmid " + plasmid + ".")

# add DeepTF lists to dataframe
df['DeepTFactor Locus_ID'] = deeptf_all
df['DeepTF Count'] = deeptf_count

# write dataframe to csv
print("Writing to csv...")
with open('plasmid_summary.csv', 'w', newline='') as f:
    df.to_csv(f, index=False)

print("Done!")

