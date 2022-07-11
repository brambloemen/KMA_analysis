.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggdendro)
library(grid)
library(gridExtra)
library(ggpubr)


source("./Import_functions.R")
source("./Plotting_functions.R")

# process zymo GMS file
zymo_GMS <- read.csv("./Zymo_GMS_D6331_plus_spike_in_I.csv", sep=";")
zymo_GMS <- zymo_GMS %>% 
  mutate(Organism = str_extract(Species,
                                "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?"),
         Total_gDNA = sum(Genomic.DNA)) %>%
  group_by(Organism) %>%
  summarize(Perc_gDNA = sum(Genomic.DNA),
            p_bpTotal = sum(Genomic.DNA)/unique(Total_gDNA),
            KMA_experiment = "zymo_GMS",
            .groups = "drop") %>%
  arrange(desc(Perc_gDNA), Organism)

# sort by descending zymo GMS gDNA% -> make organisms into ordered factors
zymo_GMS$Organism <- factor(zymo_GMS$Organism, levels = zymo_GMS$Organism,
                            ordered = TRUE)

# import and summarize files
kmafiles <- list.files() %>% unlist()
kmafiles <- str_subset(kmafiles, "(\\.mapstat|\\.res)")
kmafiles <- unique(str_remove(kmafiles, "(\\.mapstat|\\.res)"))
kmafiles <- as.list(kmafiles)

system.time(kma_all <- do.call("rbind_clean_kmaexp", kmafiles))

# append the reference community as if it were another experiment
kma_sum <- rbindlist(list(kma_all, zymo_GMS), fill = TRUE) %>%
  select(-Perc_gDNA)

# merge expected abundance onto each organism for each experiment
kma_sum <- merge(kma_sum, zymo_GMS[,c("Organism","Perc_gDNA")], by="Organism", all=TRUE)
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


# plot heatmap with dendogram of KMAexperiment
# Lots of fine-tuning required to align dendogram nicely with the heatmap
kma_heatmap <- kma_sum %>%
  filter(Organism %in% zymo_GMS$Organism,
         KMA_experiment != "zymo_GMS") %>%
  # limit anything >2x expected relative abundance
  mutate(KMA_experiment = str_remove_all(KMA_experiment, "(P|p)ure_GMS_?"),
         KMA_experiment = str_remove_all(KMA_experiment, "spike[in]?I_"),
         KMA_experiment = str_remove_all(KMA_experiment, "_cleaned"),
         KMA_experiment = str_remove_all(KMA_experiment, "(MG_|BBl_)"),
         norm_abundance= case_when(
           100*p_bpTotal/Perc_gDNA >= 2 ~ 2,
           100*p_bpTotal/Perc_gDNA <= 0 ~ 0,
           TRUE ~ 100*p_bpTotal/Perc_gDNA
         ))

# make matrix to cluster experiments
cl <- kma_heatmap %>%
  select(Organism, norm_abundance, KMA_experiment) %>%
  pivot_wider(names_from = KMA_experiment, values_from = norm_abundance, values_fill = 0)
clmatrix <- as.matrix(cl[2:ncol(cl)])
rownames(clmatrix) <- cl$Organism
clmatrix <- t(clmatrix)

# hierarchical clustering of experiments
ord <- hclust(dist(clmatrix))$order
ord_exper <- rownames(clmatrix)[ord]
# order the experiments
kma_heatmap$KMA_experiment <- ordered(kma_heatmap$KMA_experiment, levels = ord_exper)

# now using ggdendro
kma_heatmap_dendro <- as.dendrogram(hclust(d = dist(x = clmatrix)))
dendro_plot <- ggdendrogram(data = kma_heatmap_dendro, rotate = FALSE, size=5) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10))
dendro_plot

heatmap <- kma_heatmap %>%
  ggplot(aes(x=KMA_experiment, y=Organism)) +
  geom_tile(aes(fill=norm_abundance)) +
  scale_fill_gradient2(name="Relative abundance normalized to\nrelative abundance in GMS",
                       low="blue", mid="#e6e6e6", high="red", midpoint = 1) +
  scale_y_discrete(limits=rev) + #organisms with high on top, low on bottom
  # scale_x_discrete(position = "top") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
heatmap

dev.new(width=16, height=10)

print(heatmap,
      vp = viewport(x = 0.5, y = 0.275, width = 1, height = 0.55))
print(dendro_plot,
      vp = viewport(x = 0.435, y = 0.76, width = 0.5, height = 0.45))