import os
import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio import SeqIO
import seaborn as sns


#%% Import and trim alignment

alignment = AlignIO.read("../align/saccostrea/saccostrea_aligned.fasta", "fasta")

dummy = []

for i in range(len(alignment)):
    dummy.append(SeqRecord(Seq("-"),
                           id = alignment[i].id,
                           name = alignment[i].name,
                           description = alignment[i].description)
                 )

alignment_trimmed = MultipleSeqAlignment(dummy)

threshold = 0

for i in range(len(alignment[1])):
      if (alignment[:,i].count('-')) <= threshold:
          alignment_trimmed = alignment_trimmed + alignment[:,i:i+1]
        
alignment_trimmed = alignment_trimmed[:,1:]


#%% Identify haplotypes

sequences = []

for i in range(len(alignment_trimmed)):
    sequences.append(str(alignment_trimmed[i,:].seq))

haplotype_seqs = list(set(sequences))

sample_list = []
for i in range(len(alignment_trimmed)):
    sample_list.append(alignment_trimmed[i].description)

haplotype_dict = {}
for i in range(len(haplotype_seqs)):
    haplotype_dict.update({haplotype_seqs[i]: i})

sample_to_hap = pd.DataFrame({"Sample" : sample_list,
                              "Haplotype" : [""] * len(sample_list)})

for i in range(len(sample_to_hap)):
    sample_to_hap.iloc[i,1] = haplotype_dict[alignment_trimmed[i].seq]
    
#%% Export haplotype dataframe

genebank_to_lineage_df = pd.read_excel('../../Chapter_1_tables.xlsx',
                                    sheet_name = 'Tabel1')
genebank_to_lineage_dict = {}

for i in range(len(genebank_to_lineage_df)):
    genebank_to_lineage_dict.update(
        {genebank_to_lineage_df.iloc[i,1] : genebank_to_lineage_df.iloc[i,0]}
        )

def rename_ref_seqs(sample_name):
    if len(sample_name.split(" "))>1:
        return(genebank_to_lineage_dict[sample_name.split(" ")[0]])
    else:
        return(sample_name)

sample_to_hap["id"] = sample_to_hap['Sample'].apply(rename_ref_seqs)

haplotype_dict_2 = dict(zip(haplotype_dict.values(), haplotype_dict.keys()))

get_seq = lambda x: haplotype_dict_2[x]
sample_to_hap["Seq"] = sample_to_hap["Haplotype"].apply(get_seq)

sample_to_hap.to_csv("../output/sample_to_haplotype.csv")


#%% Export haplotype FASTA

haplotype_fasta = []

output_handle = open("../align/saccostrea/haplotypes.fasta", "w")

for haplotype in haplotype_dict_2:
    if " " in list(sample_to_hap['id'][sample_to_hap['Haplotype'] == haplotype])[0]:
        output_handle.write(">" + list(sample_to_hap['id'][sample_to_hap['Haplotype'] == haplotype])[0].replace(" ", "_") + "\n")
    else:
        output_handle.write(">Haplotype_" + str(haplotype) + "\n")
    output_handle.write(haplotype_dict_2[haplotype] + "\n")
    
output_handle.close()


#%% test

" " in list(sample_to_hap['id'][sample_to_hap['Haplotype'] == 3])[0]
