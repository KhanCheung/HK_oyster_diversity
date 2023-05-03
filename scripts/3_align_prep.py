import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline

#%% Get species list

species_list = pd.read_csv("../output/species_id.csv", index_col=0)

crassostrea_list = species_list[species_list['Genus'] == "Crassostrea"]
saccostrea_list = species_list[species_list['Genus'] == "Saccostrea"]

#%% Get fasta

crassostrea_records = []
saccostrea_records = []

with open("../output/oysters.fasta") as handle:
    for record in SeqIO.FastaIO.FastaIterator(handle):
        
        if record.id in crassostrea_list['Sample'].to_list():
            crassostrea_records.append(record)
        if record.id in saccostrea_list['Sample'].to_list():
            saccostrea_records.append(record)


#%% Add reference sequences to Saccostrea fasta

for seq in os.listdir("../align/saccostrea/refs"):
    
    path = os.path.join("../align/saccostrea/refs", seq)
   
    with open(path) as handle:
        for record in SeqIO.FastaIO.FastaIterator(handle):
            saccostrea_records.append(record)

SeqIO.write(saccostrea_records, "../align/saccostrea/saccostrea.fasta", "fasta")

                                