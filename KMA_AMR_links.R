.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)
library(RColorBrewer)

source("./Import_functions.R")
Refseq_results <- read.csv("./data/ResFinder_results_tab.txt", sep="\t") %>%
  mutate(Species = clean_GMS_refseq_temps(Contig),
         Resistance.gene = clean_ARG_names(Resistance.gene)) %>%
  select(-Accession.no., -Contig, -Position.in.contig) %>%
  group_by(Species, Resistance.gene) %>%
  select(names(.[sapply(., is.numeric)])) %>%
  summarise_all(sum)

Refseq_resfinder <- import_kmaexp("./data/GMSspikeI_refseq_ResFdb")
GMSspikeI_AMR <- read.csv("./data/GMSspikeI_AMRprofile.csv", sep = ";") %>%
  mutate(True_AMR=TRUE)

# process zymo GMS file
zymo_GMS <- read.csv("./data/Zymo_GMS_D6331_plus_spike_in_I.csv", sep=";")
zymo_GMS <- zymo_GMS %>% 
  mutate(Ref_species = str_extract(Species,
                                "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?"),
         Total_gDNA = sum(Genomic.DNA)) %>%
  group_by(Ref_species) %>%
  summarize(Perc_gDNA = sum(Genomic.DNA),
            p_bpTotal = sum(Genomic.DNA)/unique(Total_gDNA),
            KMA_experiment = "zymo_GMS",
            .groups = "drop") %>%
  arrange(desc(Perc_gDNA), Ref_species)

# sort by descending zymo GMS gDNA% -> make organisms into ordered factors
zymo_GMS$Ref_species <- factor(zymo_GMS$Ref_species, levels = zymo_GMS$Ref_species,
                            ordered = TRUE)



AMRlinkfiles <- list.files("./data")
AMRlinkfiles <- AMRlinkfiles[str_detect(AMRlinkfiles,"AMRlinks.csv")]
experiments <- lapply(AMRlinkfiles, str_remove, "_AMRlinks.csv")
AMRlinkfiles <- lapply(1:length(AMRlinkfiles), function(k) paste0("./data/", AMRlinkfiles[k]))

for(file in AMRlinkfiles){
  AMRdata <- Combine_AMR_data(file, "DTUdb")
  AMRdata <- merge(AMRdata, GMSspikeI_AMR, 
                   by.x=c("Ref_species", "AMR_gene") , by.y=c("Ref_species", "Resistance.gene"),
                   all = TRUE) %>%
    mutate(True_AMR=ifelse(is.na(True_AMR), FALSE, TRUE),
           AMR_coverage_depth=n_match_bases_AMR/AMR_ref_length,
           Temp_coverage_depth = bpTotal/mean_template_length)
  AMRdata <- merge(AMRdata, zymo_GMS[,c("Ref_species","Perc_gDNA")], by="Ref_species", all=TRUE)
  
  plot <- AMRdata %>%
    ggplot(aes(x = mean_templateID, y= AMR_coverage_depth, col=True_AMR)) +
    geom_point() +
    geom_text(aes(label=AMR_gene), nudge_y=-0.1) +
    scale_y_log10() +
    scale_x_log10() +
    ggtitle(str_remove_all(file, "(./data/|BBl_|_AMRlinks.csv)"))
  print(plot)
}

