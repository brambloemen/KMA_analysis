##################################################################################
#                                                                                #
#             Functions for importing and processing KMA results                 #
#                                                                                #
#                  Author: Bram Bloemen                                          #
#                  Date:   08/04/2022                                            #
#                  Date last update: 11/04/2022                                  #
#                                                                                #
##################################################################################

.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)

# function returning NCBI accession nrs (for mapping against DTU db), and organism names
# fixme: regex != precise, messy template names occur in DTU db; 
#   - maybe match against db of all species?
#   - use acession numbers to lookup organism with rentrez if stringr returns NA
# fixme: instead of organism, return more rich taxonomy data
get_orgn <- function(dataframe){
  dataframe <- dataframe %>%
    mutate(Accession = str_extract(dataframe[[1]], "[:alpha:]+_?[:alpha:]*[:digit:]+?\\.?[:digit:]+"),
           Organism=str_remove_all(dataframe[[1]], "Synthetic|\\[|\\]|\\scf.")) %>%
    mutate(Organism = str_extract(Organism, 
                                  "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?"),
           Organism = str_replace(Organism, "Clostridium difficil+e", "Clostridioides difficile")) # C. diff. was renamed
  dataframe <- as.data.frame(dataframe) #data.table sometimes behaves different to data.frame
  return(dataframe)
}

# for single experiment, format and combine the KMA .res and .mapstat files
import_kmaexp <- function(kmafile){
  kma_mapstat <- fread(file=paste0(kmafile, ".mapstat"), sep="\t", skip=6,
                       integer64 = "numeric") # otherwise numbers larger than 2^31 will lead to ggplot errors
  kma_res <- fread(file=paste0(kmafile, ".res"), sep="\t")
  
  kma <- merge(kma_mapstat, kma_res, by.x="# refSequence", by.y="#Template")
  
  kma <- get_orgn(kma)
  kma <- mutate(kma, KMA_experiment = kmafile) # add variable describing data origin
  return(kma)
}

# function to combine KMA datasets into one dataframe, and summarize data by species
rbind_clean_kmaexp <- function(..., QueryID = 80){
  kmafiles <- list(...)
  
  kma_all <- data.frame()
  for (file in kmafiles){
    kma <- import_kmaexp(file)
    # warning: this summary removes all other variables
    kma <- kma %>%
      # Calculate totals: excludes most of the unmapped reads as KMA doesn't fully output these
      mutate(Total_bp = sum(bpTotal, na.rm = TRUE),
             Total_readCount = sum(readCount, na.rm = TRUE)) %>%
      # QueryID cutoff: only templates that sufficiently matched the aligned reads are included
      filter(Query_Identity > QueryID) %>%
      group_by(Organism, KMA_experiment) %>%
      summarize(p_bpTotal = sum(bpTotal, na.rm = TRUE)/unique(Total_bp),
                mean_queryID = weighted.mean(Query_Identity, bpTotal, na.rm=TRUE),
                mean_query_coverage = weighted.mean(Query_Coverage, bpTotal, na.rm=TRUE),
                mean_templateID = weighted.mean(Template_Identity, Template_length, na.rm=TRUE),
                template_length = sum(Template_length, na.rm=TRUE),
                mean_template_coverage = weighted.mean(Template_Coverage, Template_length, na.rm=TRUE),
                bpTotal = sum(bpTotal, na.rm = TRUE),
                p_readCount = sum(readCount, na.rm = TRUE)/unique(Total_readCount),
                readCount = sum(readCount, na.rm=TRUE),
                mean_readlength = bpTotal/readCount,
                refConsensusSum = sum(refConsensusSum),
                .groups = "drop")
    
    kma_all <- rbind(kma_all, kma)
  }
  
  return(kma_all)
}



clean_GMS_refseq_temps <- function(template){
  template <- str_replace_all(template, "_", " ")
  template <- str_replace_all(template, "E.coli.*", "Escherichia coli")
  template <- str_replace_all(template, "\\..*", "")
  template <- paste0(str_to_upper(str_sub(template, 1, 1)),
                           str_sub(template, 2))
  template <- str_extract(template, "^[:alnum:]+\\s[:alnum:]+")
  template <- str_replace(template, "Clostridium difficille", "Clostridioides difficile")
  return(template)
}

clean_ARG_names <- function(amr_gene){
  amr_gene <- str_remove(amr_gene, "(_|-).+")
  return(amr_gene)
}