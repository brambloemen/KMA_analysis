.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)
library(RColorBrewer)

source("./Import_functions.R")
Refseq_results <- read.csv("ResFinder_results_tab.txt", sep="\t") %>%
  mutate(Species = clean_GMS_refseq_temps(Contig),
         Resistance.gene = clean_ARG_names(Resistance.gene)) %>%
  select(-Accession.no., -Contig, -Position.in.contig) %>%
  group_by(Species, Resistance.gene) %>%
  select(names(.[sapply(., is.numeric)])) %>%
  summarise_all(sum)


test <- fread("./BBl_pure_GMSspikeI_CB_RAD_AMRlinks.csv", fill = TRUE) %>%
  mutate(AMR_gene = clean_ARG_names(AMR_gene),
         Ref_species = get_orgn(Ref_species)) %>%
  group_by(AMR_gene, Ref_species) %>%
  summarize_all(sum)

test <- merge(test, Refseq_results, 
              by.x = c("Ref_species", "AMR_gene"), by.y=c("Species", "Resistance.gene"),
              all = TRUE)

Pure_QD_LSK <- fread("./BBl_pure_GMSspikeI_QuickDNA_LSK_AMRlinks.csv", fill = TRUE) %>%
  mutate(AMR_gene = clean_ARG_names(AMR_gene)) %>%
  group_by(AMR_gene, Ref_species) %>%
  summarize_all(sum)
Pure_QD_LSK <- merge(Pure_QD_LSK, Refseq_results, 
                 by.x = c("Ref_species", "AMR_gene"), by.y=c("Species", "Resistance.gene"),
                 all = TRUE)