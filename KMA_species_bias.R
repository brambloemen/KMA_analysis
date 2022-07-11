.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)
library(RColorBrewer)

source("./Import_functions.R")
source("./Plotting_functions.R")

# process zymo GMS file
zymo_GMS <- read.csv("./Zymo_GMS_D6331_plus_spike_in_I.csv", sep=";")
zymo_GMS <- zymo_GMS %>% 
  mutate(Organism = str_extract(Species,
                                "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?"),
         Total_gDNA = sum(Genomic.DNA)) %>%
  group_by(Organism) %>%
  mutate(Perc_gDNA = sum(Genomic.DNA),
            p_bpTotal = sum(Genomic.DNA)/unique(Total_gDNA),
            KMA_experiment = "zymo_GMS",
            .groups = "drop") %>%
  arrange(desc(Perc_gDNA), Organism)
zymo_GMS <- rbind(zymo_GMS[1:3,], zymo_GMS[8:nrow(zymo_GMS),])

# sort by descending zymo GMS gDNA% -> make organisms into ordered factors
zymo_GMS$Organism <- factor(zymo_GMS$Organism, levels = zymo_GMS$Organism,
                            ordered = TRUE)

# import and summarize files
kmafiles <- list.files() %>% unlist()
kmafiles <- str_subset(kmafiles, "(\\.mapstat|\\.res)")
kmafiles <- unique(str_remove(kmafiles, "(\\.mapstat|\\.res)"))
kmafiles <- as.list(kmafiles)

system.time(kma_all <- do.call("rbind_clean_kmaexp", kmafiles))

kma_sum <- rbindlist(list(kma_all), fill = TRUE) 
kma_sum <- merge(kma_sum, zymo_GMS[,c("Organism","Perc_gDNA","Gram.Stain")], by="Organism", all=TRUE)
kma_sum <- kma_sum %>% 
  mutate(
    Matrix=case_when(str_detect(KMA_experiment, "(P|p)ure") ~ "Pure",
                     str_detect(KMA_experiment, "(F|f)ecal") ~ "Fecal",
                     TRUE ~"zymo_GMS"),
    Kit = case_when(str_detect(KMA_experiment, "CB") ~ "CB",
                    str_detect(KMA_experiment, "QuickDNA") ~ "QuickDNA",
                    TRUE ~"zymo_GMS"),
    Library = case_when(str_detect(KMA_experiment, "RAD") ~ "RAD",
                        str_detect(KMA_experiment, "LSK") ~ "LSK",
                        TRUE ~"zymo_GMS"),
    Database = case_when(str_detect(KMA_experiment, "DTUdb") ~ "DTUdb",
                         str_detect(KMA_experiment, "GMSspikeIdb") ~ "GMSspikeIdb",
                         TRUE ~"zymo_GMS"),
    Experiment = case_when(Matrix == "zymo_GMS" ~ "zymo_GMS",
                           TRUE ~ paste(Matrix, Kit, Library, sep="_")
    ))
kma_sum <- kma_sum %>% mutate(ExperimentID = case_when(Experiment=="zymo_GMS" ~ "zymo_GMS",
                                                       TRUE ~ paste(Experiment, Database, sep="_")))


# summary by gram stain
kma_gram <- kma_sum %>%
  filter(Database == "DTUdb") %>%
  group_by(Experiment, Matrix, Gram.Stain) %>%
  summarize(p_bpTotal = sum(p_bpTotal))
zymo_gram <- zymo_GMS %>%
  group_by(Gram.Stain) %>%
  summarize(Perc_gDNA = sum(Perc_gDNA),
            p_bpTotal = Perc_gDNA/100,
            Experiment = "zymo_GMS")
kma_gram <- merge(kma_gram, zymo_gram[,c("Perc_gDNA","Gram.Stain")], by="Gram.Stain", all=TRUE) 
kma_gram <- rbindlist(list(kma_gram, zymo_gram), fill=TRUE)



png("Gramstain_exp_vs_obs.png", width=1200, height = 800, units = "px")
compare_gram <- kma_gram %>% 
  mutate(Experiment = str_remove(Experiment, "(Fecal|Pure)_"),
         Matrix = ordered(Matrix, levels=c("Pure", "Fecal"))) %>%
  filter(Experiment != "zymo_GMS" & !is.na(Gram.Stain)) %>%
  bar_exp_v_obs(Gram.Stain, 100*p_bpTotal/Perc_gDNA, Experiment) +
  ylim(NA, 2) +
  ylab("Observed / Expexted relative abundance") +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(facets = vars(Matrix))

scalefactor <- max(1200, 800)/1200
compare_gram <- compare_gram +
  theme(axis.text.x = element_text(size=20*scalefactor, angle = 0, hjust=0),
        axis.text.y = element_text(size=20*scalefactor),
        axis.title.x = element_text(size=24*scalefactor),
        axis.title.y = element_text(size=20*scalefactor),
        strip.text.x = element_text(size=24*scalefactor),
        title = element_text(size=24*scalefactor),
        legend.key = element_rect(size=1*scalefactor, color=NA),
        legend.text = element_text(size=20*scalefactor))
print(compare_gram)
dev.off()



# read-level statistics

readstats <- fread(file = "./BBl_pure_GMSspikeI_QuickDNA_LSK_DTUdb_readstats.csv") %>%
  mutate(KMA_experiment="Enzymatic - Ligation")
readstats2 <- fread(file = "./BBl_pure_GMSspikeI_QuickDNA_RAD_DTUdb_readstats.csv") %>%
  mutate(KMA_experiment="Enzymatic - Rapid")
readstats3 <- fread(file = "./BBl_pure_GMSspikeI_CB_RAD_DTUdb_readstats.csv") %>%
  mutate(KMA_experiment="Bead-beating - Rapid")
readstats <- rbind(readstats, readstats2, readstats3) %>%
  mutate(Species = ordered(Species, levels=levels(zymo_GMS$Organism)),
         KMA_experiment = ordered(KMA_experiment, levels=c("Enzymatic - Ligation", "Enzymatic - Rapid", "Bead-beating - Rapid")))
rm(readstats2, readstats3)

readstats_gram <- merge(readstats, zymo_GMS[,c("Species", "Gram.Stain")], by ="Species")
readstats_gram <- readstats_gram %>% filter(Gram.Stain != "+/-") %>%
  mutate(Gram.Stain = ifelse(Gram.Stain == "n/a", "Yeast", Gram.Stain))
  
readlength_by_spec_pure <- readstats %>%
  filter(Species %in% zymo_GMS$Organism) %>%
  ggplot(aes(x=Species, y=Length, fill=KMA_experiment)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(NA, 15000) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45))
png("readlength_by_spec_pure.png", width=1200, height = 800, units = "px")

scalefactor <- max(1200, 800)/1200
readlength_by_spec_pure <- readlength_by_spec_pure +
  theme(axis.text.x = element_text(size=20*scalefactor, angle = -45, hjust=0),
        axis.text.y = element_text(size=20*scalefactor),
        axis.title.x = element_text(size=24*scalefactor),
        axis.title.y = element_text(size=24*scalefactor),
        strip.text.x = element_text(size=24*scalefactor),
        title = element_text(size=24*scalefactor),
        legend.key = element_rect(size=1*scalefactor, color=NA),
        legend.text = element_text(size=20*scalefactor))
print(readlength_by_spec_pure)
dev.off()

mean_readlengths_by_spec <- readstats %>%
  filter(Species %in% zymo_GMS$Organism) %>%
  # mutate(KMA_experiment = as.character(KMA_experiment)) %>%
  group_by(Species, KMA_experiment) %>%
  summarize(mean_readlength = mean(Length)) %>%
  ggplot(aes(x=Species, y=mean_readlength, fill=KMA_experiment)) +
  geom_col(position = "dodge", color="black") +
  ylim(NA, 15000) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_discrete(name=NULL, labels=c("Enzymatic - Ligation", "Enzymatic - Rapid", "Bead-beating - Rapid"))
png("mean_readlengths_by_spec.png", width=1200, height = 800, units = "px")

scalefactor <- max(1200, 800)/1200
mean_readlengths_by_spec <- mean_readlengths_by_spec +
  theme(axis.text.x = element_text(size=20*scalefactor, angle = -45, hjust=0),
        axis.text.y = element_text(size=20*scalefactor),
        axis.title.x = element_text(size=24*scalefactor),
        axis.title.y = element_text(size=24*scalefactor),
        strip.text.x = element_text(size=24*scalefactor),
        title = element_text(size=24*scalefactor),
        legend.key = element_rect(size=1*scalefactor, color=NA),
        legend.text = element_text(size=20*scalefactor))
print(mean_readlengths_by_spec)
dev.off()

reatstats_by_gram_pure <- readstats_gram %>%
  ggplot(aes(x=Gram.Stain, y=Length, fill=KMA_experiment)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylim(NA, 15000) +
  xlab("Gram stain") +
  scale_fill_discrete(name=NULL)
png("reatstats_by_gram_pure.png", width=1200, height = 800, units = "px")

scalefactor <- max(1200, 800)/1200
reatstats_by_gram_pure <- reatstats_by_gram_pure +
  theme(axis.text.x = element_text(size=20*scalefactor),
        axis.text.y = element_text(size=20*scalefactor),
        axis.title.x = element_text(size=24*scalefactor),
        axis.title.y = element_text(size=24*scalefactor),
        strip.text.x = element_text(size=24*scalefactor),
        title = element_text(size=24*scalefactor),
        legend.key = element_rect(size=1*scalefactor, color=NA),
        legend.text = element_text(size=20*scalefactor))
print(reatstats_by_gram_pure)
dev.off()


# read length distribution by species
for(spec in zymo_GMS$Organism){
  filename <- paste0("violin_readlength_", spec, ".png")
  filename <- str_replace_all(filename, " ", "_")
  png(filename, width=1200, height = 800, units = "px")
  
  scalefactor <- max(1200, 800)/1200
  violin <- readstats %>%
    filter(Species == spec) %>%
    ggplot(aes(x=KMA_experiment, y=Length, fill=KMA_experiment)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45))
  violin <- violin +
    theme(axis.text.x = element_text(size=20*scalefactor, angle = -45, hjust=0),
          axis.text.y = element_text(size=20*scalefactor),
          axis.title.x = element_text(size=24*scalefactor),
          axis.title.y = element_text(size=24*scalefactor),
          strip.text.x = element_text(size=24*scalefactor),
          title = element_text(size=24*scalefactor),
          legend.key = element_rect(size=1*scalefactor, color=NA),
          legend.text = element_text(size=20*scalefactor))
  print(violin)
  dev.off()
}



# species/genus-level statistics, including unmapped reads
files <- list.files() %>% unlist()
files <- str_subset(files, "(mapstat\\.csv)")
files <- as.list(files)
mapstats <- list()
i <- 1
for(file in files){
  maps <- fread(file) %>%
    mutate(KMA_experiment=str_remove(paste(file), "_mapstat\\.csv"))
  mapstats[[i]] <- maps
  i <- i + 1
}
mapstats <- rbindlist(mapstats)

zymo_genus <- str_extract(zymo_GMS$Species, "\\w*")
mapstats <- mapstats %>% 
  mutate(
    Zymo_species = case_when(
      Species == "Unmapped" ~ "Unmapped reads",
      !(Species %in% c(as.character(zymo_GMS$Organism), "Veillonella sp")) ~ "Non-mock species",
      Species %in% c(as.character(zymo_GMS$Organism), "Veillonella sp") ~ "Mock species"),
    Zymo_species = ordered(Zymo_species,levels=c("Unmapped reads", "Non-mock species", "Mock species")),
    Genus = str_extract(Species, "\\w*"),
    Zymo_genus = case_when(
      Genus == "Unmapped" ~ "Unmapped reads",
      !(Genus %in% zymo_genus) ~ "Non-mock species",
      Genus %in% zymo_genus ~ "Mock species"),
    Zymo_genus = ordered(Zymo_genus,levels=c("Unmapped reads", "Non-mock species", "Mock species"))) %>% 
  mutate(
    Matrix=case_when(str_detect(KMA_experiment, "(P|p)ure") ~ "Pure Mock",
                     str_detect(KMA_experiment, "(F|f)ecal") ~ "Mock spiked in Fecal",
                     TRUE ~"zymo_GMS"),
    Kit = case_when(str_detect(KMA_experiment, "CB") ~ "Bead-beating",
                    str_detect(KMA_experiment, "QuickDNA") ~ "Enzymatic",
                    TRUE ~"zymo_GMS"),
    Library = case_when(str_detect(KMA_experiment, "RAD") ~ "Rapid",
                        str_detect(KMA_experiment, "LSK") ~ "Ligation",
                        TRUE ~"zymo_GMS"),
    Database = case_when(str_detect(KMA_experiment, "DTUdb") ~ "DTU database",
                         str_detect(KMA_experiment, "GMSspikeIdb") ~ "Mock community database",
                         TRUE ~"zymo_GMS"),
    Experiment = case_when(Matrix == "zymo_GMS" ~ "zymo_GMS",
                           TRUE ~ paste(Matrix, Kit, Library, sep="_")),
    Matrix = ordered(Matrix, levels=c("Pure Mock", "Mock spiked in Fecal")),
    Database = ordered(Database, levels=c("Mock community database", "DTU database", "zymo_GMS")))

mapstats %>% 
  group_by(KMA_experiment, Matrix, Database) %>%
  mutate(Total=sum(Total_readlength),
         Method = paste0(Kit, " - ", Library), 
         Method = ordered(Method, levels = c("Enzymatic - Ligation", "Enzymatic - Rapid", "Bead-beating - Rapid"))) %>%
  group_by(Method, Matrix, Database, Zymo_species) %>%
  summarize(Total_readlength = sum(Total_readlength)/unique(Total)) %>%
  ggplot() +
  geom_col(aes(x=Method, y=Total_readlength, fill=Zymo_species), position = "stack") +
  scale_y_continuous(name="Percent of bases sequenced", labels = scales::percent) +
  theme(axis.text.x = element_text(angle=-45, hjust=0),
        legend.position = "top") + 
  scale_fill_brewer(palette="Set1", name=NULL) +
  facet_grid(rows = vars(Matrix), cols= vars(Database))

mapstats %>% 
  group_by(KMA_experiment, Matrix, Database) %>%
  mutate(Total=sum(Total_readlength),
         Method = paste0(Kit, " - ", Library), 
         Method = ordered(Method, levels = c("Enzymatic - Ligation", "Enzymatic - Rapid", "Bead-beating - Rapid"))) %>%
  group_by(Method, Matrix, Database, Zymo_species) %>%
  summarize(Total_readlength = sum(Total_readlength)) %>%
  ggplot() +
  geom_col(aes(x=Method, y=as.numeric(Total_readlength), fill=Zymo_species), position = "stack") +
  ylab("Bases sequenced") +
  theme(axis.text.x = element_text(angle=-45, hjust=0),
        legend.position = "top") +
  scale_fill_brewer(palette="Set1", name=NULL) +
  facet_grid(rows = vars(Matrix), cols= vars(Database))

NonGMS <- mapstats %>% filter(Matrix == "GMS spiked in Fecal", Database=="DTU database", Zymo_species == "Non-GMS species")
