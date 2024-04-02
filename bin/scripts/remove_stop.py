# remove_stop.py
# Remove stop characters '*' from FASTA file
# Input: FASTA file
# Output: FASTA file without stop characters



file_in = snakemake.input[0]
file_out = snakemake.output[0]

# Read in FASTA file
print("Reading in file %s" % file_in)
with open(file_in,'r') as f:
    data = f.read()
    # Get length of file read
    print("Length of file read: %s" % len(data))

# Remove stop characters '*'
print("Removing stop characters '*'")
data = data.replace('*','')
# Get length of file after removing '*'
print("Length of file after removing '*': %s" % len(data))

# Write to output file
print("Writing file to %s" % file_out)
with open(file_out,'w') as f:
    f.write(data)
