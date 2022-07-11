.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)
library(caret)

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
zymo_genus <- str_extract(zymo_GMS$Organism, "\\w*")
zymo_organisms <- zymo_GMS$Organism


Pure_QD_LSK <- import_kmaexp("BBl_pure_GMSspikeI_QuickDNA_LSK_DTUdb") %>%
  get_orgn() %>%
  mutate(Genus = str_extract(Organism, "\\w*"))
Pure_QD_LSK <- merge(Pure_QD_LSK, zymo_GMS[,c("Organism","Perc_gDNA")], by="Organism", all=TRUE) %>%
  mutate(is_GMS_genus = ifelse(Genus %in% zymo_genus, TRUE, FALSE),
         is_GMS_organism = ifelse(Organism %in% zymo_organisms, TRUE, FALSE))%>%
  # Calculate totals: excludes most of the unmapped reads as KMA doesn't fully output these
  mutate(Total_bp = sum(bpTotal, na.rm = TRUE),
         p_bpTotal = bpTotal/Total_bp)

Pure_QD_LSK <- Pure_QD_LSK %>% filter(Query_Identity>0)
mod_o <- glm(is_GMS_organism ~ log10(p_bpTotal) + Template_Coverage + Query_Identity + Template_Identity + Depth + Template_length, 
           data = Pure_QD_LSK, family=binomial)
Pure_QD_LSK$Prediction_o <- predict(mod_o, Pure_QD_LSK, type="response")
ggplot(Pure_QD_LSK, aes(x=p_bpTotal, y=Prediction_o, color=is_GMS_organism)) + 
  geom_point() + 
  scale_x_log10(limits=c(1e-8,0.3), labels = scales::percent)
ggplot(Pure_QD_LSK, aes(x=p_bpTotal, y=Prediction_o, color=is_GMS_genus)) + 
  geom_point() + 
  scale_x_log10(limits=c(1e-8,0.3), labels = scales::percent)



Pure_CB_RAD <- import_kmaexp("BBl_pure_GMSspikeI_CB_RAD_DTUdb") %>%
  get_orgn() %>%
  mutate(Genus = str_extract(Organism, "\\w*"))
Pure_CB_RAD <- merge(Pure_CB_RAD, zymo_GMS[,c("Organism","Perc_gDNA")], by="Organism", all=TRUE) %>%
  mutate(is_GMS_genus = ifelse(Genus %in% zymo_genus, TRUE, FALSE),
         is_GMS_organism = ifelse(Organism %in% zymo_organisms, TRUE, FALSE))%>%
  # Calculate totals: excludes most of the unmapped reads as KMA doesn't fully output these
  mutate(Total_bp = sum(bpTotal, na.rm = TRUE),
         p_bpTotal = bpTotal/Total_bp)
Pure_QD_RAD <- import_kmaexp("BBl_pure_GMSspikeI_QuickDNA_RAD_DTUdb") %>%
  get_orgn() %>%
  mutate(Genus = str_extract(Organism, "\\w*"))
Pure_QD_RAD <- merge(Pure_QD_RAD, zymo_GMS[,c("Organism","Perc_gDNA")], by="Organism", all=TRUE) %>%
  mutate(is_GMS_genus = ifelse(Genus %in% zymo_genus, TRUE, FALSE),
         is_GMS_organism = ifelse(Organism %in% zymo_organisms, TRUE, FALSE))%>%
  # Calculate totals: excludes most of the unmapped reads as KMA doesn't fully output these
  mutate(Total_bp = sum(bpTotal, na.rm = TRUE),
         p_bpTotal = bpTotal/Total_bp) %>%
  mutate(MG_criteria = ordered(case_when(
    Depth >= 20 & Template_Coverage >= 75 & Template_length >=1e05 & Template_Identity >=90 ~ "8 High +",
    Depth >= 15 & Template_Coverage >= 75 & Template_length >=1e05 & Template_Identity >=80 ~ "7 High -",
    Depth >= 10 & Template_Coverage >= 75 & Template_length >=5e04 & Template_Identity >=80 ~ "6 Normal +",
    Depth >= 5 & Template_Coverage >= 75 & Template_length >=1e04 & Template_Identity >=67.5 ~ "5 Normal -",
    Depth >= 2 & Template_Coverage >= 50 & Template_length >=1e03 & Template_Identity >=60 ~ "4 Low +",
    Depth >= 1 & Template_Coverage >= 50 & Template_length >=1e03 & Template_Identity >=45 ~ "3 Low -",
    Depth >= 0.5 & Template_Identity >=40 ~ "2 Very Low",
    Template_Identity >= 0 ~ "1 Trace",
    TRUE ~ "0 Absence"
  )))


Pure_QD_RAD <- Pure_QD_RAD %>% filter(Query_Identity>0)
Pure_QD_RAD$Prediction_o <- predict(mod_o, Pure_QD_RAD, type="response")
Pure_QD_RAD <- Pure_QD_RAD %>% 
  mutate(Predict_GMS_organism = Prediction_o > 0.5,
         Correct_prediction_o = Predict_GMS_organism == is_GMS_organism,
         Correct_prediction_g = Predict_GMS_organism == is_GMS_genus)
ggplot(Pure_QD_RAD, aes(x=p_bpTotal, y=Prediction_o, color=is_GMS_organism)) + 
  geom_point() + 
  scale_x_log10(limits=c(1e-8,0.3), labels = scales::percent)
ggplot(Pure_QD_RAD, aes(x=p_bpTotal, y=Prediction_o, color=is_GMS_genus)) + 
  geom_point() + 
  scale_x_log10(limits=c(1e-8,0.3), labels = scales::percent)
ggplot(Pure_QD_RAD, aes(x=p_bpTotal, y=Prediction_o, color=MG_criteria)) + 
  geom_point() + 
  scale_x_log10(limits=c(1e-8,0.3), labels = scales::percent) +
  scale_color_brewer(palette = "RdYlGn")




# import and summarize files
kmafiles <- list.files() %>% unlist()
kmafiles <- str_subset(kmafiles, "(\\.mapstat|\\.res)")
kmafiles <- unique(str_remove(kmafiles, "(\\.mapstat|\\.res)"))
kmafiles <- kmafiles[kmafiles != "GMSspikeI_refseq_ResFdb"]
kmafiles <- kmafiles[str_detect(kmafiles, "(P|p)ure")]
kmafiles <- kmafiles[str_detect(kmafiles, "DTUdb")]
kmafiles <- as.list(kmafiles)

system.time(kma_all <- do.call("rbind_clean_kmaexp", kmafiles))


# merge expected abundance onto each organism for each experiment
kma_sum <- merge(kma_all, zymo_GMS[,c("Organism","Perc_gDNA")], by="Organism", all=TRUE)
kma_sum <- kma_sum %>% mutate(is_GMS_species = ifelse(!is.na(Perc_gDNA), TRUE, FALSE))