#! python3

import os
import pandas as pd

contig_list=[]

def get_contig(tab):
    df=pd.read_csv(tab,sep='\t',usecols=['#FILE','SEQUENCE']) #opens tab file as pd.df
    _n=df.shape[0] #get number of rows in df
    n=0
    
    while n<_n: #will only run for number of rows in pd.df
        file=str(df['#FILE'][n]) #looks at reach row 'n' of pd.df
        contig='>'+str(df['SEQUENCE'][n]) #looks at corresponding contig code
        
        if contig in contig_list: #if contig was already added to list, continues to next contig
            n=n+1
            continue
        else:
            contig_list.append(contig)

        with open(file,'r') as f: #opens file identified in pd.df
            data=f.read()
            start=data.find(contig) #starts at contig identified in pd.df
            end=data.find('>', start+1) #determines end of contig with '>'
            result=data[start:end] #isolates contig from '>' from '>'

            contig_directory='contigs_plasmid'
            if not os.path.isdir(contig_directory):
                os.mkdir(contig_directory)

            _fname=str(df['#FILE'][n])[:-15]
            fname=_fname+'_'+str(df['SEQUENCE'][n])+'.fasta' #final file name example = 0flye_contig_1.fasta
            path=contig_directory+'/'+fname #contigs/0flye_contig_1.fasta
            print(path)
            with open(snakemake.output[0], 'a') as f:
                f.write(result) #outputs selected contigs to appropriate .fasta file
            f.close()

        n=n+1

if __name__=="__main__":
    get_contig(snakemake.input[0])