import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline 
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

#%% Make Blast database

cline = NcbimakeblastdbCommandline(dbtype = "nucl",
                                   input_file = "../blast/blastdb/ostreidae.fasta",
                                   out = "../blast/blastdb/ostreidae")

cline()


#%% Blastn

cline = NcbiblastnCommandline(query = "../output/oysters.fasta",
                              db = "../blast/blastdb/ostreidae",
                              evalue = 0.001, outfmt=5, max_target_seqs = 1,
                              out = "../blast/blast_out.xml",
                              num_threads = 8)

cline()


#%% Parse Blast output

blast_records = NCBIXML.parse(open("../blast/blast_out.xml"))

sample_list = []
binomial_list = []


for blast_record in blast_records:
    
    sample = blast_record.query

    if blast_record.alignments:
        title = blast_record.alignments[0].title
        binomial = title.split("|")[3]
    else:
        binomial = np.NaN
        
    sample_list.append(sample)
    binomial_list.append(binomial)
    

#%% Construct species ID DataFrame

df = pd.DataFrame(
    {'Sample': sample_list,
     'Binomial': binomial_list
     })
 
df.dropna(subset=['Binomial'], inplace=True)
df.reset_index(inplace=True)

site_list = []
genus_list = []
species_list = []

def modify_entry(entry):
    site_list.append(entry['Sample'].split("_")[0])
    
    if (len(entry['Binomial'].split(" ")) == 1):
        genus_list.append(entry['Binomial'])
        species_list.append("sp.")
    else:
        genus_list.append(entry['Binomial'].split(" ")[0])
        species_list.append(entry['Binomial'].split(" ")[1])
    
df.apply(modify_entry, axis=1)


df.insert(loc=0, column='Site', value = pd.Series(site_list))
df.insert(loc=3, column='Genus', value = pd.Series(genus_list))
df.insert(loc=4, column='Species', value = pd.Series(species_list))

df.drop(["index", "Binomial"], axis = 1, inplace = True)


#%% Export species ID dataframe as csv

df.to_csv("../output/species_id.csv")
