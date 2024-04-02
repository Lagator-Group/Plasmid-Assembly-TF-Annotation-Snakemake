# Extract contigs with plasmid markers from a tab-separated file with columns #FILE and SEQUENCE
# and write to a new fasta file with name and sequence.

import os
import pandas as pd

contig_list = []  # list to store contig codes (e.g. '>Contig1') to prevent duplicate contigs being extracted


def get_contig(tab):
    """Extract contigs from a tab-separated file with columns #FILE and SEQUENCE
    and write to a new fasta file with name and sequence.

    Parameters
    ----------
    tab : str
        Path to tab-separated file.
    """
    print(f"Reading tab file: {tab}")
    df = pd.read_csv(tab, sep='\t', usecols=['#FILE', 'SEQUENCE'])
    print(f"Number of rows in df: {df.shape[0]}")
    _n = df.shape[0]  # number of rows in df
    n = 0  # counter for while loop
    while n < _n:  # will only run for number of rows in pd.df
        print(f"Processing row {n+1} of {_n}")
        file = str(df['#FILE'][n])  # filename from pd.df
        contig = '>' + str(df['SEQUENCE'][n])  # contig code from pd.df
        if contig in contig_list:  # if contig was already added to list, continue to next contig
            print(f"Contig '{contig}' already extracted, skipping...")
            n = n + 1
            continue
        else:
            contig_list.append(contig)
            print(f"Adding contig '{contig}' to contig_list")
        with open(file, 'r') as f:  # open file identified in pd.df
            data = f.read()
            start = data.find(contig)  # start of contig in file
            end = data.find('>', start+1)  # end of contig with '>'
            result = data[start:end]  # isolate contig from '>' to '>'
            contig_directory = 'contigs_plasmid'
            _fname = str(df['#FILE'][n])[:-15]  # filename without extension
            fname = _fname + '_' + str(df['SEQUENCE'][n]) + '.fasta'
            path = contig_directory + '/' + fname
            with open(snakemake.output[0], 'a') as f:  # write to output file
                print(f"Writing contig to {path}")
                f.write(result)
            f.close()
        n = n + 1


if __name__ == "__main__":
    get_contig(snakemake.input[0])  # run get_contig function on input file



