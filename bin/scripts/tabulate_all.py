# %%
import pandas as pd
import os

# %%
plasmid_list = []
for plasmid in os.listdir('fasta_plasmid'):
    plasmid = plasmid.replace('.fasta','')
    plasmid_list.append(plasmid)

# %%
df=pd.DataFrame(columns=['Plasmid'])
df['Plasmid'] = plasmid_list

# %%
plasmid_length = []
for plasmid in df['Plasmid']:
    with open('fasta_plasmid/'+plasmid+'.fasta') as f:
        data = f.readlines()
        length = len(data[1])
    plasmid_length.append(length)

df['Plasmid Length'] = plasmid_length

# %%
resistance_all = []
gene_all = []
resistance_count = []
for plasmid in df['Plasmid']:
    print(plasmid)
    resistance_list = []
    gene_list = []
    amr_df = pd.read_csv('abricate_amr/'+plasmid+'.tab',sep='\t')
    if len(amr_df) == 0:
        gene_list.append('None')
        resistance_list.append('None')
        resistance_count.append(0)
    else:
        for gene in amr_df['GENE']:
            gene_list.append(gene)
        for resistance in amr_df['RESISTANCE']:
            resistance_list.append(resistance)
        resistance_count.append(len(resistance_list))
    resistance_all.append(resistance_list)
    gene_all.append(gene_list)

df['Resistance'] = resistance_all
df['Resistance Gene'] = gene_all
df['Resistance Count'] = resistance_count

# %%
tf_prot_all = []
tf_gene_all = []
tf_count = []
for plasmid in df['Plasmid']:
    print(plasmid)
    tf_prot_list = []
    tf_gene_list = []
    tf_df = pd.read_csv('sprot_TF/'+plasmid+'.tsv',sep='\t')
    if len(tf_df) == 0:
        tf_gene_list.append('None')
        tf_prot_list.append('None')
        tf_count.append(0)
    else:
        for gene in tf_df['Gene Names']:
            tf_gene_list.append(gene)
        for prot in tf_df['Protein names']:
            tf_prot_list.append(prot)
        tf_count.append(len(tf_prot_list))
    tf_prot_all.append(tf_prot_list)
    tf_gene_all.append(tf_gene_list)

tf_gene_simplified = []
for gene_list in tf_gene_all:
    simplified_list = []
    for gene in gene_list:
        gene=str(gene)
        simplified_list.append(gene.split(' ')[0])
    tf_gene_simplified.append(simplified_list)

df['TF Protein'] = tf_prot_all
df['TF Gene'] = tf_gene_all
df['TF Gene Simplified'] = tf_gene_simplified
df['TF Count'] = tf_count

# %%
deeptf_all = []
deeptf_count = []
for plasmid in df['Plasmid']:
    print(plasmid)
    deeptf_list = []
    deeptf_df = pd.read_csv('deepTFactor/'+plasmid+'/prediction_result.txt',sep='\t')
    n=0
    for prediction in deeptf_df['prediction']:
        prediction = str(prediction)
        if prediction == 'True':
            locid = deeptf_df['sequence_ID'][n]
            deeptf_list.append(locid)
        else:
            pass
        n+=1
    deeptf_count.append(len(deeptf_list))
    deeptf_all.append(deeptf_list)

df['DeepTFactor Locus_ID'] = deeptf_all
df['DeepTF Count'] = deeptf_count

# %%
with open('plasmid_summary.csv', 'w', newline='') as f:
    df.to_csv(f, index=False)