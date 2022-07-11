#############################################################
#                                                           #
#                  Analysis of KMA results                  #
#                                                           #
#                  Author: Bram Bloemen                     #
#                  Date:   21/02/2022                       #
#                  Date last update: 26/02/2022             #
#                                                           #
#############################################################
# requirements:
#   - csv file of Zymo D6331 gut microbiome standard in working dir: Zymo_GMS_D6331_composition.csv
#   - KMA output files: .res and .mapstat (produced by kma option -ef), located in working dir

# outputs:
#   - histogram of read counts vs expected for pure GMS
#   - histogram of total bases sequenced vs expected for pure GMS

# future improvements:
#   - improve species parsing; more descriptive taxonomy (Veillonella rogosae: not present in DTU db?)
#   - include log plot for very low abundance species
#   - adapt for use on linux server -> where to install R packages?
#   - make callable from command line/bash script

.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
# library(BiocManager)
# library(rentrez)
library(data.table)
# library(zeallot)

zymo_GMS <- read.csv("./Zymo_GMS_D6331_composition.csv", sep=";")
zymo_GMS <- mutate(zymo_GMS, Organism_wo_strain = str_extract(Species, 
                                         "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?"))
# because multiple Ecoli strains are present in zymo GMS:
zymo_GMS_unique <- zymo_GMS %>% group_by(Organism_wo_strain) %>%
  summarize(Perc_gDNA = sum(Genomic.DNA))


# function to process the necessary data from mapstat and res files
# produces histograms
GMS_species_barplots <- function(kmafile, QueryID=80){
  
  kma_mapstat_fp <- paste0(kmafile, ".mapstat")
  kma_res_fp <- paste0(kmafile, ".res")
  kma_mapstat <- fread(file=kma_mapstat_fp, sep="\t", skip=6)
  kma_res <- fread(file=kma_res_fp, sep="\t")
  
  # function using regex to return NCBI accession nrs and organism names from mapstat
  # fixme: regex != precise; 
  #   - maybe match against db of all species?
  #   - use acession numbers to lookup organism with rentrez if stringr returns NA
  # fixme: instead of organism, return entire taxonomic tree
  get_orgn <- function(dataframe){
    
    dataframe <- dataframe %>%
      mutate(Accession = str_extract(dataframe[[1]], "^([:alnum:]|\\.)*"),
             Organism=str_remove(dataframe[[1]], "Synthetic")) %>%
      mutate(Organism = str_extract(Organism, 
                                    "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?"))
    
    return(dataframe)
  }
  
  # combine mapstat and res file
  kma <- get_orgn(kma_mapstat)
  kma <- merge(kma, kma_res, by.x="# refSequence", by.y="#Template")
  
  # filter out Query ID < QueryID
  kma <- filter(kma, Query_Identity > QueryID)
  
  # summarize and prepare dataframe for plotting
  kma_sum <- 
    merge(kma, zymo_GMS_unique,
          by.x="Organism", by.y="Organism_wo_strain", all=TRUE) %>%
    mutate(Total_reads=sum(readCount, na.rm = TRUE),
           Total_bp = sum(bpTotal, na.rm = TRUE)) %>% 
    group_by(Organism) %>%
    summarize(readCount = sum(readCount, na.rm = TRUE)/unique(Total_reads), #/unique(), to only keep 1 value per group
              bpTotal = sum(bpTotal, na.rm = TRUE)/unique(Total_bp),
              Perc_gDNA = unique(Perc_gDNA)/100)

  # filter out species from zymo GMS only
  # Fixme: create another "species" with all alignments to non-GMS species
  kma_sum <- filter(kma_sum, Organism %in% zymo_GMS$Organism_wo_strain)
  
  # sort by descending zymo GMS gDNA% -> make organisms into ordered factors
  kma_sum <- kma_sum[order(-kma_sum$Perc_gDNA, kma_sum$Organism), ]
  kma_sum$Organism <- factor(kma_sum$Organism, levels = kma_sum$Organism,
                             ordered = TRUE)
  
  # dataframe for plotting
  kma_sum_long <- pivot_longer(data= kma_sum, 
                               cols=2:4,
                               names_to = "Variable",
                               values_to = "Value") %>%
    mutate(Variable = factor(Variable, c("readCount", "bpTotal", "Perc_gDNA"),
                             ordered = TRUE))
  kma_sum_long <- kma_sum_long[order(kma_sum_long$Organism),]
  
  # function to plot histograms
  # Fixme: indicate low abundant species with * mark on x-axis
  plot_hists <- function(dataframe, var_to_plot){
    hist <- dataframe %>%
      filter(Variable %in% c(var_to_plot, "Perc_gDNA")) %>%
      ggplot(aes(x=Organism, y=Value, group=Variable)) + 
      geom_col(position="dodge", aes(fill=Variable)) +
      theme_bw() +
      scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_text(angle=90),
            panel.grid.major.x = element_blank(),
            legend.title = element_blank()) +
      scale_fill_manual(values=c("#00BFC4", "#F8766D"), labels = c("Observed", "Expected for\npure GMS"))
    return(hist)
  }
  
  hist_bpTotal <- plot_hists(kma_sum_long, "bpTotal") +
    labs(x = "Species", y = "Percentage of bases")
  return(hist_bpTotal)
}

# get filenames for working directory
kmafiles <- list.files() %>% unlist()
kmafiles <- str_subset(kmafiles, "(.*mapstat|.*res)")
kmafiles <- unique(str_remove(kmafiles, "(\\.mapstat|\\.res)"))

# loop over files, make histograms
# fixme: make new plots for low abundance species: use the data included in ggplot list objects
# fixme: replace for loop with parallel lapply
for (file in kmafiles){
  hist_bpTotal <- GMS_species_barplots(file)
  
  hist_bpTotal <- hist_bpTotal +
    labs(title = "Distribution of mapped bases") +
    scale_y_continuous(limits = c(NA, 0.25), labels = scales::percent)
  png(filename=paste0(file, "_hist_bpTotal.png"), width = 1000, height=500)
  print(hist_bpTotal)
  dev.off()
}
