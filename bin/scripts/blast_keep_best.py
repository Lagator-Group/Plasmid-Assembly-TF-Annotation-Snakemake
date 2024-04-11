# Script to extract the best hit from a BLASTx output file
# for each locus tag in the file
# Output is a tab-delimited file with two columns:
#    Locus Tag and Entry
# If no blast hits are found for a locus tag,
# output is a single-row DataFrame with a message

import pandas as pd

# Read the .tsv file into a DataFrame
# Read the .tsv file into a DataFrame
print("Reading input file...")
try:
    # Read the blastx output file into a Pandas DataFrame
    blastx_df = pd.read_csv(snakemake.input[0],
                            header=None,  # no header in the file
                            delimiter='\t')  # tab-delimited file
    print("Input file read successfully.")

    # Extract the locus tags and entries from the DataFrame
    locus_list = blastx_df[0].to_list()  # locus tag column
    entry_list = blastx_df[1].to_list()  # entry name column
    print("Extracted locus and entry lists from input file.")

    # Create lists to store the best locus and entry for each locus tag
    best_locus = []
    best_entry = []
    print("Created empty lists for best locus and entry.")

    # Iterate over the locus tags in the file
    n = 0  # track the row number
    for locus in locus_list:
        print("Processing locus {}...".format(locus))
        # If this locus tag is not already in the best lists,
        # add it and its corresponding entry to the lists
        if locus in best_locus:
            pass
        else:
            best_locus.append(locus)
            best_entry.append(entry_list[n])
            print("Added best hit for locus {}: {}".format(locus, entry_list[n]))

        # Increment the row number
        n += 1

    # Create a new DataFrame with the best locus and entry for each locus tag
    best_df = pd.DataFrame(
        {'Locus Tag': best_locus,  # column names
         'Entry': best_entry})
    print("Created final DataFrame with best hits.")

except:
    # If there was an error in reading the file, create a DataFrame with a single row
    best_df = pd.DataFrame({'Locus Tag': ['No BLAST results'],  # column names
                            'Entry': ['No BLAST results']})
    print("Error occurred while reading input file:", sys.exc_info()[0])

# Write the best hits to a new file
with open(snakemake.output[0], 'w', newline='') as f:
    best_df.to_csv(f, index=False, sep='\t')  # write the DataFrame to file
    print("Wrote best hits to output file.")

