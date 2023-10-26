# %%
import pandas as pd

# %%
trembl_df=pd.read_csv(snakemake.input[0],sep='\t', header=None)
column_headers = ['Locus Tag','Entry','2','3','4','5','6','7','8','9','10','11']
trembl_df.columns = column_headers

# %%
with open(snakemake.output[0],'w') as f:
    trembl_df.to_csv(f,index=False,sep='\t')


