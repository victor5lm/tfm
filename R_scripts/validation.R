# This R script describes how validation of the best model obtained from
# the HEI-based classification of the study cohort was performed
# using the relative abundances table from the validation cohort.

library(pROC)
library(phyloseq)
library(tidyr)
library(dplyr)
library(mikropml)
library(viridis)
library(RColorBrewer)
library(readr)

# Auxiliary function for ROC curve plotting:

roc.plot <- function(roc.list, title.in, p.ci = FALSE,
                     shuffle = TRUE){
  # generates palette
  if (shuffle == TRUE){
    lengthroc <- length(roc.list)/2
    pal <- rep(viridis::viridis(lengthroc), 2)
    aes_ggroc <- c("linetype", "colour")
    line_types <- c(rep("solid", lengthroc), rep("twodash", lengthroc))
  }
  else {
    pal <- viridis::viridis(length(roc.list)) # Graph colour can be changed here
    aes_ggroc <- "colour"
    line_types <- "solid"
  }
  
  # extracts auc
  roc.list %>% 
    map(~tibble(AUC = .x$auc)) %>% 
    bind_rows(.id = "name") -> data.auc
  
  # generates labels
  data.auc %>% 
    mutate(label_long=paste0(name,", AUC = ",paste(round(AUC,1))),
           label_AUC=paste0("AUC = ",paste(round(AUC,2)))) -> data.labels
  
  names(roc.list) <- data.labels$label_long
  
  p <- ggroc(c(roc.list),
             size = 1.2, legacy.axes = TRUE, aes = aes_ggroc) + 
    geom_line(size = 1.2) +
    labs(x = "False Positive Rate", 
         y = "True Positive Rate",
         title = title.in) +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = line_types) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted") +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "bottom") +
    coord_equal()
  
  if (p.ci == TRUE){
    ## plots + confidence intervals
    ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, l = 25))
    
    dat.ci.list <- lapply(ci.list, function(ciobj) 
      data.frame(x = as.numeric(rownames(ciobj)),
                 upper =  1- ciobj[, 1],
                 lower = 1- ciobj[, 3]))
    
    for(i in 1:length(roc.list)) {
      p <- p + geom_ribbon(
        data = 1 - dat.ci.list[[i]],
        aes(x = x, ymin = lower, ymax = upper),
        fill = pal[i],
        alpha = 0.2,
        inherit.aes = FALSE) 
    }}
  return(p)
}

# First, we need the training data used by the best model
# trained on the relative abundances table from the 
# study cohort grouped on their HEI.
training_data_rf_mikropml <- originals_hei_rf$"88"$trained_model$trainingData

# Now, we need the relative abundances table from the validation cohort:
abundances_val_clr_hei <- read.csv("clr_abundances_validation_cohort.csv")

# Once we have done this, we will use this abundances table for validation
# asking the ML classifier to predict classes on this new data.
# To do so, however, the abundances table from the validation cohort
# needs to have the same columns as the training data used by 
# our model. 

# We will create a new variable containing all columns shared by 
# both dataframes and the relative clr-transformed abundances
# from the validation cohort.
common_col_names <- insersect(names(training_data_rf_mikropml), names(abundances_val_clr_hei))

# Next, from the abundances table from the validation cohort, 
# we will keep only these common columns
df1 <- abundances_val_clr_hei[, common_col_names, drop = FALSE]

# Next, we need to find those columns (or taxa) from the training data that are
# not in the abundances table from the validation cohort
diff_col_names <- setdiff(names(training_data_rf_mikropml), names(abundances_val_clr_hei))

# We then create a new dataframe containing all taxa from the study cohort
# not shared with the validation cohort
df2 <- training_data_rf_mikropml[, diff_col_names, drop = FALSE]

# We then need to change the rownames
rownames(df2) <- NULL
rownames(df2) <- seq_len(nrow(df2))

# Those taxa that do not appear on the abundances table from the
# validation cohort will have the lowest CLR-abundance
# registered in this table, indicating that these taxa
# were not abundant in the samples from this cohort
df2[] <- -3.005435

# We make both dataframes have the same number of rows
df2 <- df2[-tail(seq_len(nrow(df2)), 3), , drop = FALSE]

# Finally, we can merge both dataframes: the one containing all taxa
# from the training data not shared with the validation cohort followed
# by the other dataframe with all taxa shared by both cohorts.
# Please remember, since the whole final dataframe will contain 
# CLR-transformed relative abundances relative to the validation cohort,
# all abundances from the first dataframe will have the lowest CLR-value,
# since they were not registered when doing the taxonomic assignment of the
# validation cohort, and all abundances from the second dataframe will contain 
# all CLR-transformed relative abundances from the validation cohort
# relative to all taxa shared by both cohorts
final_df <- cbind(df2, df1)

# We then add a final column indicating the HEI group each sample from the
# validation cohort belongs to
final_df$HEI_group <- abundances_val_clr_hei$HEI_group

# We then calculate the area under the curve
roc_validation <- roc(as.factor(final_df$HEI_group), (predict(originals_ias_rf$`88`$trained_model, final_df, type = "prob"))$Good_HEI)

# And finally plot the ROC curve relative to the performance of our best model
# on this validation dataset
my.roc_val <- list("Validation cohort" = roc_validation)
plot_roc_validation <- roc.plot(my.roc_val, "", p.ci = TRUE, shuffle = FALSE)
