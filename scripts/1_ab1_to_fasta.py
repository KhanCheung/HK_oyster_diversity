import glob
import numpy as np
import pandas as pd
from Bio import SeqIO

raw_dir = "../raw/" 

output_dir = "../output/"

    
#%% Generate dataframe for sample metadata

ab1_paths = glob.glob(raw_dir + '/*/*.ab1', recursive=True)
ab1_filename = list(i.split("\\")[-1] for i in ab1_paths)
site = [i.split("_")[0] for i in ab1_filename]
replicate = [i.split("_")[1] for i in ab1_filename]
sample = list(i+"_"+j for i,j in zip(site, replicate))
replicate = [int(x) for x in replicate]

sample_list = pd.DataFrame(
    {"Sample": sample,
     "Site": site,
     "Replicate": replicate,
     "ab1_path": ab1_paths
    }
)


#%% Parse ab1 files to fasta format

fasta_records = []

def abi_to_fasta(entry):
        
    input_handle = open(entry["ab1_path"], "rb")
    
    record = next(SeqIO.AbiIO.AbiIterator(input_handle, trim=True))
    
    if (len(record.seq) >= 600):
        record.id = record.name.split("__")[0]
        fasta_records.append(record)

    input_handle.close()
    
    
sample_list.apply(abi_to_fasta, axis=1)

SeqIO.write(fasta_records, "../output/oysters.fasta", "fasta")
