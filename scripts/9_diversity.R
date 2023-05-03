library("rstudioapi")
library('dplyr')
library("iNEXT")
library('vegan')

setwd(dirname(getActiveDocumentContext()$path))


# Import and rearrange data ----------------------------------------------------

raw_data = read.csv('../output/species_id_updated.csv')[ , -1]

for (i in 1:nrow(raw_data)){
  if (raw_data[i,3] == "Crassostrea"){
    raw_data[i,3] = "Magallana"
  }
}

raw_data['Binomial'] = paste(raw_data$Genus, raw_data$Species)

species = as.list(unique(raw_data['Binomial']))
sites = as.list(unique(raw_data['Site']))

species_comp = matrix(0, nrow = lengths(species), ncol = lengths(sites))

colnames(species_comp) = unlist(sites)
rownames(species_comp) = sort(unlist(species))

for (site in colnames(species_comp)){
  for (species in rownames(species_comp)){
    
    species_comp[species, site] = dim(raw_data[raw_data$Site == site & raw_data$Binomial == species,])[1]
    
  }
}

write.csv(species_comp, "../output/species_matrix.csv")

# Sample coverage esit --------------------------------------------------------------

diversity_est = iNEXT(species_comp, q = c(0,1,2))

ggiNEXT(diversity_est, type=2, facet.var = "Assemblage")


# Beta diversity ---------------------------------------------------------------

chisq.test(species_comp, simulate.p.value = TRUE)

species_dist = vegdist(t(species_comp))
mds <- metaMDS(species_dist)
mds_data <- as.data.frame(mds$points)
mds_data$SampleID <- rownames(mds_data)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = SampleID)) +
  geom_point()

