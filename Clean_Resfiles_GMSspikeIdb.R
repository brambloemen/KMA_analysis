# cleans data of KMA with GMS as reference database
# GMS reference database contains non-scientific species IDs in fasta headers
#   => also in KMA result files
# --> convert unclear naming to scientific names

.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)

# input: common name of KMA output files
Clean_KMA_mapto_GMS <- function(filename){
  
  # Read in mapstat and res files
  df_mapstat <- fread(file=paste0(filename, ".mapstat"), sep="\t", skip=6)
  df_mapstat_header <- fread(file=paste0(filename, ".mapstat"), sep="\t", nrows=6, header = FALSE)
  df_res <- fread(file=paste0(filename, ".res"), sep="\t")
  
  replace_templates <- function(dataframe){
    dataframe[[1]] <- str_replace_all(dataframe[[1]], "_", " ")
    dataframe[[1]] <- str_replace_all(dataframe[[1]], "E.coli.*", "Escherichia coli")
    dataframe[[1]] <- str_replace_all(dataframe[[1]], "\\..*", "")
    dataframe[[1]] <- paste0(str_to_upper(str_sub(dataframe[[1]], 1, 1)),
                             str_sub(dataframe[[1]], 2))
    dataframe[[1]] <- str_extract(dataframe[[1]], "^[:alnum:]+\\s[:alnum:]+")
    dataframe[[1]] <- str_replace(dataframe[[1]], "Clostridium difficille", "Clostridioides difficile")
    # summarise results per species:
    # merging mapstat to res file in KMA_analysis.R requires only 1 row per template species
    #   however, multiple Ecoli strains are present --> summarise to Ecoli on species level
    #   cur_data_all() gives the current data for the current group (incl. grouping vars)
    #   --> used for referring to first column, regardless of column/variable names
    dataframe <- dataframe %>% group_by(select(cur_data_all(),1)) %>%
      summarise_all(sum)
    return(dataframe)
  }
  
  df_mapstat <- replace_templates(df_mapstat)
  write_tsv(df_mapstat_header, paste0(filename,"_cleaned.mapstat"))
  write_tsv(df_mapstat, paste0(filename,"_cleaned.mapstat"), append = TRUE, col_names = TRUE)
  
  df_res <- replace_templates(df_res)
  write_tsv(df_res, paste0(filename,"_cleaned.res"))
}

files <- list.files() 
files <- files[str_detect(files,"GMSspikeIdb(\\.mapstat|\\.res)")]
files <- files[str_detect(files,"cleaned(\\.mapstat|\\.res)", negate = TRUE)]
files <- str_remove(files, "\\.mapstat|\\.res")
files <- unique(files)
files <- as.list(files)

lapply(files, Clean_KMA_mapto_GMS)
# Clean_KMA_mapto_GMS(files[[1]])
# Clean_KMA_mapto_GMS(files[[2]])
