# Scripts for the generation of the most accurate ML classifier for the HPF classification of the Study cohort:

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(phyloseq)
library(viridis)
library(RColorBrewer)
library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(purrr)
library(mikropml)

# Classifier construction

library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(viridis)

if( !require("mikropml") ){
    install.packages("mikropml")
    library(mikropml)
}

data_complete <- read.table('data_clr_complete.tsv',
                      header = TRUE, sep = '\t',# skip = 1,
                      row.names = 1)
ntrees<-seq(1,100)
tuning_rf<-list(mtry=seq(1,70))
originals <- list()

for (ntree in ntrees){print(ntree)
    set.seed(2019)
    results <- run_ml(data_complete, "rf",
    outcome_colname = 'Grupo_HPF',
    cross_val=caret::trainControl(method = "LOOCV"),
    training_frac = 0.80, seed=2019,
    calculate_performance = TRUE, ntree=ntree,
    hyperparameters = tuning_rf,
    find_feature_importance = TRUE)
                      
    key<-toString(ntree)
    originals[[key]]<-results
}

# Best model
originals$`17`$performance$AUC #0.88333

# Feature importances
imp.bi <- varImp(originals$`17`$trained_model, scale = FALSE)

# Heatmap for the 20 most important microbial taxa identified by the best ML model trained after the HPF classification of the study cohort
a <- imp.bi$importance %>% 
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    arrange(desc(Overall)) %>% mutate(rowname = forcats::fct_inorder(rowname)) %>%
    rename(feature = rowname) %>%
    mutate(type = "binary")
a_imp <- a[1:20,]

groups_imp <- c("0","0","0","1","1","0","0","0","0","0","1","1","0","0","0","1","0","1","1","1")
groups_imp <- as.data.frame(groups_imp)
a_imp <- cbind(a_imp, groups_imp$groups_imp)

a_imp_ordered <- a_imp[order(a_imp$"groups_imp$groups_imp"),]
groups_for_pheatmap_genuses <- read.csv("heatmap_files/additional_files_hpf_classification/groups_for_pheatmap_genuses.csv")
groups_for_pheatmap <- read.csv("heatmap_files/additional_files_hpf_classification/groups_for_pheatmap.csv")
abunda_clr_pheatmap <- read.csv("heatmap_files/additional_files_hpf_classification/abundancias_clr_para_pheatmap.csv")
rownames(abunda_clr_pheatmap) <- abunda_clr_pheatmap$Sample
abunda_clr_pheatmap$Sample <- NULL
abunda_clr_pheatmap_matrix <- as.matrix(abunda_clr_pheatmap)

p2 = ComplexHeatmap::HeatmapAnnotation(Importance = anno_barplot(a_imp_ordered$Overall, gp = gpar(fill=ifelse(a_imp_ordered$`groups_imp$groups_imp`==0,"#073444","#bae9fa"))), height = unit(3, "cm"), HPF_Group=groups_for_pheatmap_genuses$Group,col=list(HPF_Group=c("Low HPF consumption (≤ 15%)"="#073444","High HPF consumption (> 15%)"="#bae9fa")),annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), annotation_legend_param = list(title="", labels_gp = gpar(font=2, fontsize = 12)), border = TRUE, gap = unit(1, "mm"))

p3 = ComplexHeatmap::rowAnnotation(HPF_Group=groups_for_pheatmap$Group,col=list(HPF_Group=c("Low HPF consumption (≤ 15%)"="#073444","High HPF consumption (> 15%)"="#bae9fa")),annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), annotation_legend_param = list(title="", labels_gp = gpar(font=2, fontsize = 12)), border = TRUE, show_legend = FALSE, show_annotation_name = FALSE)

heatmap_pheatmat <- ComplexHeatmap::Heatmap(abunda_clr_pheatmap_matrix, cluster_rows = FALSE, top_annotation = p2, right_annotation = p3, col = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n = 9)), cluster_columns = FALSE, show_row_names = FALSE, heatmap_legend_param = list(title="", legend_height = unit(4,"cm"), labels_gp = gpar(font=2), legend_direction = "vertical"), column_names_gp = grid::gpar(fontsize = 10.5, fontface = "bold.italic"), rect_gp = gpar(col = "#073444", lwd = 1), row_split = groups_for_pheatmap$Group, row_gap = unit(2, "mm"), row_title_gp = gpar(col = c("white","white")), column_split = a_imp_ordered$`groups_imp$groups_imp`, column_title_gp = gpar(col = c("white","white")), column_gap = unit(2, "mm"))

draw(heatmap_pheatmat,heatmap_legend_side = "left", annotation_legend_side = "bottom")
