##################################################################################
#                                                                                #
#                  Comparison of KMA results between experiments                 #
#                                                                                #
#                  Author: Bram Bloemen                                          #
#                  Date:   26/02/2022                                            #
#                  Date last update: 11/07/2022                                  #
#                                                                                #
##################################################################################

.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)
library(RColorBrewer)


source("./Import_functions.R")
source("./Plotting_functions.R")


# process zymo GMS file
zymo_GMS <- read.csv("./data/Zymo_GMS_D6331_plus_spike_in_I.csv", sep=";")
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
kmafiles <- list.files("./data") %>% unlist()
kmafiles <- str_subset(kmafiles, "(\\.mapstat|\\.res)")
kmafiles <- unique(str_remove(kmafiles, "(\\.mapstat|\\.res)"))
kmafiles <- kmafiles[kmafiles != "GMSspikeI_refseq_ResFdb"]
kmafiles <- as.list(kmafiles)
kmafiles <- lapply(1:length(kmafiles), function(k) paste0("./data/", kmafiles[[k]]))

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
    Kit = case_when(str_detect(KMA_experiment, "CB_beads") ~ "bead-beating, beads",
                    str_detect(KMA_experiment, "CB") ~ "bead-beating, column",
                    str_detect(KMA_experiment, "QuickDNA") ~ "Enzymatic, beads",
                    TRUE ~"zymo_GMS"),
    Library = case_when(str_detect(KMA_experiment, "RAD") ~ "Rapid",
                        str_detect(KMA_experiment, "LSK") ~ "Ligation",
                        TRUE ~"zymo_GMS"),
    Database = case_when(str_detect(KMA_experiment, "DTUdb") ~ "DTUdb",
                         str_detect(KMA_experiment, "GMSspikeIdb") ~ "GMSspikeIdb",
                         TRUE ~"zymo_GMS"),
    Experiment = case_when(Matrix == "zymo_GMS" ~ "zymo_GMS",
                           TRUE ~ paste(Matrix, Kit, Library, sep="_")),
    ExperimentID = case_when(Experiment=="zymo_GMS" ~ "zymo_GMS",
                             TRUE ~ paste(Experiment, Database, sep="_")))


zymo_GMS_lowab <- zymo_GMS %>%
  filter(Perc_gDNA < 0.5)
plot <- kma_sum %>% 
  filter(Perc_gDNA < 0.5) %>%
  filter(Database == "DTUdb" | Database == "zymo_GMS") %>%
  scatter_exp_v_obs(p_bpTotal, Matrix, "zymo_GMS") + 
  ylab("Observed % of mapped bases") +
  xlab("theoretical distribution") +
  theme(legend.position = "right") +
  scale_y_log10(limits=c(1e-7,NA)) + # , labels = scales::percent
  scale_x_log10(limits=c(1e-7,NA)) +
  geom_text(data=zymo_GMS_lowab,
            aes(Perc_gDNA/100, Perc_gDNA/100, label=str_replace(Organism, " ", "\n")),
            size=6,
            position = position_nudge(y=0.4, x=-0.3))
export_png("output/Low_ab_DTUdb.png", plot, width = 1200, height=800)


# plots expected vs observed, pure vs fecal per workflow
patterns <- c("(GMSspikeIdb|CB|RAD)", "(GMSspikeIdb|CB|LSK)", "(GMSspikeIdb|QuickDNA|LSK|beads)", "(GMSspikeIdb|QuickDNA|LSK|CB_RAD)")
names <- c("scatter_DTUdb_QDLSK.png", "scatter_DTUdb_QDRAD.png", "scatter_DTUdb_CBRAD.png", "scatter_DTUdb_CBbeads.png")
titles <- c("Enzymatic (zymo QuickDNA) - Ligation sequencing", "Enzymatic (zymo QuickDNA) - Rapid sequencing", 
            "Bead beating (Claremont Bio) - Rapid sequencing", "Bead beating (Claremont Bio + zymo beads) - Rapid sequencing")
for(i in 1:length(patterns)){
  scatter <- kma_sum %>% 
    filter(KMA_experiment %in% kmafiles[str_detect(kmafiles, patterns[i], negate = TRUE)] | KMA_experiment == "zymo_GMS") %>%
    scatter_exp_v_obs(p_bpTotal, Matrix, "zymo_GMS") +
    scale_y_log10(limits=c(0.000001,0.3), labels = scales::percent) +
    scale_x_log10(limits=c(0.000001,0.3), labels = scales::percent) +
    ylab("Observed % of mapped bases") +
    xlab("theoretical distribution") +
    stat_smooth(kma_sum[kma_sum$KMA_experiment %in% kmafiles[str_detect(kmafiles, patterns[i], negate = TRUE)]], 
                mapping=aes(Perc_gDNA/100, p_bpTotal, color=Matrix, fill=Matrix), 
                method = "lm", alpha=0.1, fullrange = TRUE) +
    ggtitle(titles[i])
  print(scatter)
  export_png(paste0("output/",names[i]), scatter, width = 1600, height = 1000)
}

# bar plots for inter-species comparison
for(matrix in c("Pure", "Fecal")){
  for(db in c("DTUdb", "GMSspikeIdb")){
    
    barplot <- paste0("output/Bar_", matrix, db, ".png")
    barplot_norm <- paste0("output/Bar_norm_", matrix, db, ".png")
    barplot_readlengths <- paste0("output/Bar_readlength_", matrix, db, ".png")
    barplot_totalbp <- paste0("output/Bar_totalbp_", matrix, db, ".png")
    
    df <- kma_sum %>% 
      filter(Matrix %in% c(matrix, "zymo_GMS"), 
             Organism %in% zymo_GMS$Organism, 
             Database %in% c(db, "zymo_GMS"))
    
    bars_GMS_pure <-  df %>% 
      bar_exp_v_obs(Organism, p_bpTotal, Experiment, "zymo_GMS") + 
      ylab("Percentage of bases classified")
    export_png(barplot, bars_GMS_pure, width=1600, height = 1200)
    
    bars_GMS_norm_pure <- df %>% 
      filter(Experiment != "zymo_GMS") %>%
      bar_exp_v_obs(Organism, (100*p_bpTotal/Perc_gDNA), Experiment) + 
      ylab("Relative abundance normalized \nto expected abundance") +
      ylim(NA, 3) +
      geom_hline(yintercept=1, linetype="dashed")
    export_png(barplot_norm, bars_GMS_norm_pure, width=1600, height = 1200)

    bar_readlengths <- df %>%
      bar_exp_v_obs(Organism, mean_readlength, Experiment, refcategory = "zymo_GMS")
    export_png(barplot_readlengths, bar_readlengths, width=1600, height = 1200)

    bar_totalbp <- df %>%
      bar_exp_v_obs(Organism, bpTotal, Experiment, refcategory = "zymo_GMS")
    export_png(barplot_totalbp, bar_totalbp, width=1600, height = 1200)

  }
}

zymo_GMS_species_percent <- zymo_GMS %>%
  arrange(desc(Perc_gDNA), Organism) %>%
  mutate(Species_percent = ordered(paste0(Organism, " ", signif(Perc_gDNA, digits = 3), " %"))) %>%
  select(Species_percent)
zymo_GMS_species_percent <- c(zymo_GMS_species_percent$Species_percent)

kma_sum <- kma_sum %>% 
  mutate(Species_percent = ordered(paste0(Organism, " ", signif(Perc_gDNA, digits = 3), " %"),
                                   levels=zymo_GMS_species_percent))
# Heatmaps
for(db in c("DTUdb", "GMSspikeIdb")){
  df <- kma_sum %>% filter(Database %in% c(db, "zymo_GMS"))
  filename <- paste0("output/Heatmap_", db, ".png")
  hm <- heatmap_exp_v_obs(df, Experiment, Species_percent, p_bpTotal, "zymo_GMS") +
    theme(axis.text.x = element_text(angle=-45, hjust=+1.0, vjust=0.25, size=20),
          legend.position = "right")
  export_png(name = filename, hm)
}

for(matrix in c("Pure", "Fecal")){
  df <- kma_sum %>% filter(Matrix %in% c(matrix, "zymo_GMS"))
  filename <- paste0("output/Heatmap_", matrix, ".png")
  hm <- heatmap_exp_v_obs(df, ExperimentID, Species_percent, p_bpTotal, "zymo_GMS") +
    theme(axis.text.x = element_text(angle=-45, hjust=+1.0, vjust=0.25, size=20),
          legend.position = "right")
  export_png(name = filename, hm)
  
  for(db in c("DTUdb", "GMSspikeIdb")){
    df2 <- df %>% filter(Database %in% c(db, "zymo_GMS"))
    filename <- paste0("output/Heatmap_", matrix, db, ".png")
    hm <- heatmap_exp_v_obs(df2, Experiment, Species_percent, p_bpTotal, "zymo_GMS") +
      theme(axis.text.x = element_text(angle=-45, hjust=+1.0, vjust=0.25, size=20),
            legend.position = "right")
    export_png(name = filename, hm)
  }
}

hm <- heatmap_exp_v_obs(kma_sum, ExperimentID, Species_percent, p_bpTotal, "zymo_GMS") +
  theme(axis.text.x = element_text(angle=-45, hjust=+1.0, vjust=0.25, size=20),
        legend.position = "right")
export_png(name = "output/Heatmap.png", hm)


Fecal_background <- rbind_clean_kmaexp("./raw_data/Fecal_background_QD_LSK_DTUdb") %>% 
  mutate(Matrix="Fecal_background", Kit = "QuickDNA", Library = "LSK", Database = "DTUdb",
         Experiment = paste(Matrix, Kit, Library, sep="_"), ExperimentID = paste(Experiment, Database, sep="_"))

kma_sum <- rbind(kma_sum, Fecal_background, fill=TRUE)
kma_sum_genus <- kma_sum %>% mutate(Genus = str_extract(Organism, "\\w*")) %>%
  filter(Database %in% c("zymo_GMS", "DTUdb")) %>%
  group_by(KMA_experiment, ExperimentID, Database, Matrix, Genus) %>%
  summarize(bpTotal = sum(bpTotal), Perc_gDNA = sum(Perc_gDNA)) %>%
  ungroup()
heatmap_total <- kma_sum_genus %>%
  heatmap_all_org(ExperimentID, Genus, bpTotal, refcat = "zymo_GMS", refvar = Perc_gDNA)
png("output/heatmap_total_genus.png", width = 8000, height = 8000, units = "px")
scalefactor <- 8000/1200
print(heatmap_total +
        theme(axis.text.x = element_text(size=16*scalefactor, hjust = 0, vjust = 0, angle = 45),
              axis.text.y = element_text(size=6*scalefactor),
              axis.title.x = element_text(size=24*scalefactor),
              # axis.title.y = element_text(size=24*scalefactor),
              strip.text.x = element_text(size=24*scalefactor),
              title = element_text(size=24*scalefactor),
              legend.key = element_rect(size=1*scalefactor, color=NA),
              legend.text = element_text(size=20*scalefactor),
              legend.key.size = unit(scalefactor,"cm"),
              plot.margin = margin(t=5, unit = "cm"))
)
dev.off()




kma_fecal <- kma_sum_genus %>% filter(KMA_experiment != "zymo_GMS" & Database != "GMSspikeIdb" & Matrix == "Fecal")
