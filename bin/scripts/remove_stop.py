file_in = snakemake.input[0]
file_out = snakemake.output[0]

with open(file_in,'r') as f:
    data = f.read()
    data = data.replace('*','')

with open(file_out,'w') as f:
    f.write(data)


