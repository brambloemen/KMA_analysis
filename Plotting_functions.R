.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)
library(RColorBrewer)

##################################################################
#                                                                #
#             Functions for plotting KMA results                 #
#                                                                #
#               Author: Bram Bloemen                             #
#               Date:   11/04/2022                               #
#               Date last update: 20/04/2022                     #
#                                                                #
##################################################################


# function to make scatterplot of observed value vs expected distribution
# to see if experimental observations deviates from expected
scatter_exp_v_obs <- function(df, y, categorical, refcategory){
  # df: data, needs to include a Perc_gDNA expected abundance for each expected species
  # x is always the same: expected relative abundance
  # y is quantative variable that shows observed distribution
  # categorical: category to compare (e.g. experiment)
  # refcategory: "Ground truth" category, against which to compare
  
  # double curly brackets: pass arguments as expressions, only works for dplyr-related
  data <- df %>% filter({{categorical}} != refcategory) 
  plot <- 
    ggplot() + 
    geom_point(data = data, 
               aes(x=Perc_gDNA/100, y={{y}}, color={{categorical}}), 
               position = position_dodge2(width = 0.2), size=4) +
    geom_abline(aes(slope=1, intercept=0), linetype=2) +
    scale_color_brewer(palette="Set1") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=270))
  plot <- plot + 
    scale_y_continuous(trans='log10', labels = scales::percent) +
    scale_x_continuous(name="Theoretical distribution",
                       trans='log10', labels = scales::percent)
  
  return(plot)
}



# barplots: To compare organism abundance between a limited amount of experiments
bar_exp_v_obs <- function(df, x, y, categorical, refcategory=NULL, colorpal="Set1"){
  # x: Organism variable. Can be at any taxonomic level, as long as it is present in reference community.
  # y: Quantification variable. E.g. read count, mapped bases, coverage etc.
  # categorical: variable to compare (e.g. experiments)
  # refdf: optional reference data
  # refcategory: optional reference category.
  
  # initial filter step + format fill variables
  df <- df %>%
    mutate(fillvar = {{categorical}},
           xvar = {{x}})
  df$fillvar <- ordered(df$fillvar)
  fill_lvls <- levels(df$fillvar)
  
  # check if reference category is present, remove from automatically colored categories
  if(length(refcategory!=0)){
    refdf <- df %>% filter({{categorical}} == refcategory)
    df <- df %>% 
      filter(xvar %in% refdf[[enexpr(x)]]) # enexpr returns the naked expression
    # always plot reference category as last, overwrite fill levels
    fill_lvls <- levels(df$fillvar)[levels(df$fillvar)!=refcategory]
    df$fillvar <- ordered(df$fillvar, levels=c(fill_lvls, refcategory))
  }
  
  # create color palette for non-reference levels
  fillcolors <- brewer.pal(n=length(fill_lvls), colorpal) 
  
  plot <- df %>%
    ggplot(aes(x={{x}}, y={{y}}, fill=fillvar)) +
    # here, the position_dodge2 makes sure that bars are of constant width, even when missing data
    geom_col(position=position_dodge2(preserve = "single"), color="black", width = 0.8) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle=-45, hjust=0.0),
          axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          plot.margin = margin(2, 5.5, 2, 2, "cm")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  
  # Test if ref cat is present
  if(length(refcategory)==0){
    plot <- plot + scale_fill_manual(values = fillcolors)
  }
  else{
    # Plot reference in grey
    plot <- plot + scale_fill_manual(values = c(fillcolors, "#e6e6e6"))
  }
  
  return(plot)
}



# heatmap clustering different experimental parameters
heatmap_exp_v_obs <- function(df, x, y, hue, refcat){
  # df: data, needs to include a Perc_gDNA var, with expected abundances
  # x: Categorical to be compared (e.g. experimental parameter)
  # y: Organism variable, is plotted on vertical axis
  # hue: Quantative variable that is used to fill the heatmap
  # normalise the RA of each species to RA in GMS
  # ref category: organisms in spike-in community and their expected abundances
  ref_org_filt <- df %>% filter({{x}} == refcat) %>%
    select({{y}}) %>% distinct()
  ref_org_filt <- as.character(ref_org_filt[[1]])
  
  df <- df %>%
    mutate(xaxis = {{x}}, yaxis={{y}}) %>% 
    filter(yaxis %in% ref_org_filt & !is.na(yaxis),
           xaxis != refcat) %>%
    # limit anything >2x expected relative abundance
    mutate(norm_abundance=log2(100*{{hue}}/Perc_gDNA))
      # norm_abundance= case_when(
      #   100*{{hue}}/Perc_gDNA >= 2 ~ 2,
      #   100*{{hue}}/Perc_gDNA <= 0 ~ 0,
      #   TRUE ~ 100*{{hue}}/Perc_gDNA
      # ))
  
  # make matrix to cluster experiments
  cl <- df %>%
    select(yaxis, norm_abundance, xaxis) %>%
    pivot_wider(names_from = xaxis, values_from = norm_abundance, values_fill = 0)
  clmatrix <- as.matrix(cl[2:ncol(cl)])
  rownames(clmatrix) <- cl$yaxis
  clmatrix <- t(clmatrix)
  
  # hierarchical clustering of experiments
  ord <- hclust(dist(clmatrix))$order
  ord_exper <- rownames(clmatrix)[ord]
  # order the experiments
  df$xaxis <- ordered(df$xaxis, levels = ord_exper)
  
  heatmap <- df %>%
    ggplot(aes(x=xaxis, y=yaxis)) +
    geom_tile(aes(fill=norm_abundance)) +
    scale_fill_gradient2(name="Log2 of ratio observed/expected",
                         low="blue", mid="#e6e6e6", high="red", midpoint = 0) +
    #reverse organism order on heatmap y axis (i.e. top=high abundance, bottom=low)s
    scale_y_discrete(limits=rev, name=NULL) + 
    scale_x_discrete(position = "top", name=NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=-45, hjust=-0.05, vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  return(heatmap)
}


# heatmap clustering different experimental parameters
heatmap_all_org <- function(df, x, y, hue, refcat=NULL, refvar=NULL){
  # df: data, needs to include a Perc_gDNA var, with expected abundances
  # x: Categorical to be compared (e.g. experimental parameter)
  # y: Organism variable, is plotted on vertical axis
  # hue: Quantative variable that is used to fill the heatmap
  # refcat: reference category if present
  # refvar: reference var by which organisms in refcat should be sorted
  
  
  # normalise the RA of each species
  df <- df %>%
    mutate(xaxis = {{x}}, yaxis={{y}}) %>%
    group_by(xaxis) %>%
    mutate(huetotal = sum({{hue}})) %>%
    # limit anything >2x expected relative abundance
    mutate(
      norm_abundance= {{hue}}/huetotal,
      norm_abundance=log10(norm_abundance))
  
  if(length(refcat != 0)){
    total_refvar <- df %>% filter(xaxis == refcat) %>% mutate(sortvar = {{refvar}}) %>% 
      summarise(sum(sortvar))
    total_refvar <- as.numeric(total_refvar)
    df <- df %>% mutate(sortvar = {{refvar}})
  }
  
  # make matrix to cluster samples (= experiments)
  cl <- df %>%
    select(yaxis, norm_abundance, xaxis) %>%
    pivot_wider(names_from = xaxis, values_from = norm_abundance, values_fill = 0)
  clmatrix <- as.matrix(cl[2:ncol(cl)])
  rownames(clmatrix) <- cl$yaxis
  clmatrix <- t(clmatrix)
  
  # hierarchical clustering of experiments
  ord <- hclust(dist(clmatrix))$order
  ord_exper <- rownames(clmatrix)[ord]
  # order the experiments
  df$xaxis <- ordered(df$xaxis, levels = ord_exper)
  
  # order Organisms by median occurence
  org_order <- df %>% group_by(yaxis) %>%
    summarize(med = median(norm_abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>% select(yaxis) %>% distinct()
  org_order <- as.character(org_order[[1]])
  
  # if reference category is present, correct the order: reference organisms on top
  if(length(refcat != 0)){
    org_order_ref <- df %>% ungroup() %>% filter(xaxis == refcat) %>% 
      arrange(desc(sortvar)) %>% select(yaxis) %>% distinct()
    org_order_ref <- as.character(org_order_ref[[1]])
    
    org_order_other <- df %>% ungroup() %>% filter(xaxis != refcat) %>% group_by(yaxis) %>%
      summarize(med = median(norm_abundance, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(med)) %>% select(yaxis) %>% distinct()
    org_order_other <- as.character(org_order_other[[1]])
    org_order_other <- org_order_other[!(org_order_other %in% org_order_ref)] # remove ref organisms out of other organisms
    
    org_order <- c(org_order_ref, org_order_other)
  }
  
  df <- df %>% mutate(yaxis= ordered(yaxis, levels=org_order))
  
  # define midpoint for heatmap color hue
  midpoint <- median(df$norm_abundance, na.rm = TRUE)
  
  heatmap <- df %>%
    ggplot(aes(x=xaxis, y=yaxis)) +
    geom_tile(aes(fill=norm_abundance)) +
    scale_fill_gradient2(name="log10 Relative abundance",
                         low="blue", mid="#e6e6e6", high="red", midpoint = midpoint) +
    #reverse organism order on heatmap y axis (i.e. top=high abundance, bottom=low)s
    scale_y_discrete(limits=rev, name=NULL) + 
    scale_x_discrete(position = "top", name=NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=-45, hjust=-0.5, vjust=-0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  return(heatmap)
}

export_png <- function(name, plot, width=1200, height=1200){
  png(name, width=width, height = height, units = "px")
  
  scalefactor <- max(width, height)/1200
  plot <- plot +
    theme(axis.text.x = element_text(size=20*scalefactor),
          axis.text.y = element_text(size=20*scalefactor),
          axis.title.x = element_text(size=24*scalefactor),
          axis.title.y = element_text(size=24*scalefactor),
          strip.text.x = element_text(size=24*scalefactor),
          title = element_text(size=24*scalefactor),
          legend.key = element_rect(size=1*scalefactor, color=NA),
          legend.text = element_text(size=20*scalefactor))
  print(plot)
  dev.off()
}
