# Gets sprot data for predicted locus tags

import pandas as pd

# Read in prediction file from deeptFactor
prediction_file = snakemake.input[0]

# Read in sprot file from snakemake
sprot_file = snakemake.input[1]

# Output file
out_file = snakemake.output[0]

# Read in deeptFactor prediction file
prediction_result_df = pd.read_csv(prediction_file, sep='\t')

# Iterate through deeptFactor prediction results
# and extract only those that are predicted
n=0
predicted_locid_list = []
for prediction in prediction_result_df['prediction']:
    
    # Convert prediction to string for comparison
    prediction = str(prediction)
    
    # If prediction is True, extract locus tag from sequence_ID column
    if prediction == 'True':
        locid = prediction_result_df['sequence_ID'][n]
        predicted_locid_list.append(locid)
        print(f'Locus tag found in deeptFactor prediction results: {locid}')
    else:
        pass
    
    n+=1

# Create lists for storing data
locid_list = []
entry_list = []
protein_name_list = []
gene_name_list = []

# Iterate through locus tags and retrieve information from sprot file
for locid in predicted_locid_list:
    
    # Print locus tag for troubleshooting
    print(f'Locus tag being searched in sprot file: {locid}')
    
    # Try to find locus tag in sprot file
    try:
        sprot_df = pd.read_csv(sprot_file, delimiter='\t')
        
        # Retrieve entry, protein name, and gene name from sprot file
        entry = sprot_df.loc[sprot_df['Locus Tag']==locid]['Entry'].iloc[0]
        protein_name = sprot_df.loc[sprot_df['Locus Tag']==locid]['Protein names'].iloc[0]
        gene_name = sprot_df.loc[sprot_df['Locus Tag']==locid]['Gene Names'].iloc[0]
        
        # Append data to lists
        locid_list.append(locid)
        entry_list.append(entry)
        protein_name_list.append(protein_name)
        gene_name_list.append(gene_name)
        
    # If locus tag not found in sprot file, append "Not in sprot database" to lists
    except:
        locid_list.append(locid)
        entry_list.append('Not in sprot database')
        protein_name_list.append('Not in sprot database')
        gene_name_list.append('Not in sprot database')
        print(f'Locus tag not found in sprot file: {locid}')

# Create a dataframe from the lists of data
predicted_df = pd.DataFrame({'Locus Tag':locid_list,'Entry':entry_list,'Protein names':protein_name_list,'Gene Names':gene_name_list})

# Write dataframe to output file
with open(out_file,'w',newline='') as f:
    predicted_df.to_csv(f,index=False,sep=',')



