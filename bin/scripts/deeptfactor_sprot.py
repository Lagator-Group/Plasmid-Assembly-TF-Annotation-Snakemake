# %%
import pandas as pd

# %%
plasmid = 'AP019704.1'
prediction_file = snakemake.input[0]
sprot_file = snakemake.input[1]
out_file = snakemake.output[0]

# %%
prediction_result_df = pd.read_csv(prediction_file, sep='\t')

# %%
n=0
predicted_locid_list = []
for prediction in prediction_result_df['prediction']:
    prediction = str(prediction)
    if prediction == 'True':
        locid = prediction_result_df['sequence_ID'][n]
        predicted_locid_list.append(locid)
    else:
        pass
    n+=1


# %%
locid_list = []
entry_list = []
protein_name_list = []
gene_name_list = []

for locid in predicted_locid_list:
    print(locid)
    try:
        sprot_df = pd.read_csv(sprot_file, delimiter='\t')
        entry = sprot_df.loc[sprot_df['Locus Tag']==locid]['Entry'].iloc[0]
        protein_name = sprot_df.loc[sprot_df['Locus Tag']==locid]['Protein names'].iloc[0]
        gene_name = sprot_df.loc[sprot_df['Locus Tag']==locid]['Gene Names'].iloc[0]

        locid_list.append(locid)
        entry_list.append(entry)
        protein_name_list.append(protein_name)
        gene_name_list.append(gene_name)
    except:
        locid_list.append(locid)
        entry_list.append('Not in sprot database')
        protein_name_list.append('Not in sprot database')
        gene_name_list.append('Not in sprot database')

predicted_df = pd.DataFrame({'Locus Tag':locid_list,'Entry':entry_list,'Protein names':protein_name_list,'Gene Names':gene_name_list})

# %%
with open(out_file,'w',newline='') as f:
    predicted_df.to_csv(f,index=False,sep=',')


