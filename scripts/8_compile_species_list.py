import pandas as pd

#%% Import data

df = pd.read_csv("../output/species_id.csv", index_col=0)

sac_sample_to_hap = pd.read_csv("../output/sample_to_haplotype.csv")
sac_hap_to_line = pd.read_csv("../output/saccostrea_haplotype_to_lineage.csv")

#%% Saccostrea sample to lineage

saccostrea_df = sac_sample_to_hap[['Sample', 'Haplotype']].join(sac_hap_to_line.set_index('Haplotype'),on="Haplotype")

saccostrea_dict = dict(zip(saccostrea_df["Sample"], saccostrea_df["Lineage"]))

#%% Update main df

for i in range(len(df)):
    if df.iloc[i,2] == "Saccostrea":
        df.iloc[i,3] = saccostrea_dict[df.iloc[i,1]]
        
#%%

df.to_csv("../output/species_id_updated.csv")
