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

#Classifier construction
library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(viridis)

if( !require("mikropml") ){
    install.packages("mikropml")
    library(mikropml)
}

data_complete_hei_groups_rf <- read.table('DATA/data_clr_complete_hei_groups.tsv',
                      header = TRUE, sep = '\t', # skip = 1,
                      row.names = 1)

descrCorr <- cor(data_complete_hei_groups_rf)
highCorr <- findCorrelation(descrCorr, 0.90)
data_rf.uncor <- data_complete_hei_groups_rf[, -highCorr]

ntrees <- seq(1, 100)
tuning_rf <- list(mtry = seq(1, 100))
originals_hei_rf <- list()

for (ntree in ntrees){print(ntree)
    set.seed(2019)
    results <- run_ml(data_rf.uncor, "rf",
    outcome_colname = 'HEI_group',
    cross_val = caret::trainControl(method = "LOOCV"),
    training_frac = 0.80, seed = 2019,
    calculate_performance = TRUE, ntree = ntree,
    hyperparameters = tuning_rf,
    find_feature_importance = FALSE)
    key <- toString(ntree)
    originals_hei_rf[[key]] <- results
}

#Best model
originals_hei_rf$"88"$performance #AUC = 0.85

#Variable importances
imp.bi <- varImp(originals_hei_rf$"88", scale = FALSE)

#Heatmap with the 20 most important features according to this ML classifier and based on the HEI classification of the Study cohort
table_imp_mikropml_rf <- read.csv("DATA/importances_HEI/imp_mikropml_rf_hei.csv")
table_imp_mikropml_rf_for_pheatmap <- cbind(table_imp_mikropml_rf,c(1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,1))
colnames(table_imp_mikropml_rf_for_pheatmap)[4] <- "groups_imp$groups_imp"
table_imp_mikropml_rf_for_pheatmap <- table_imp_mikropml_rf_for_pheatmap[order(table_imp_mikropml_rf_for_pheatmap$`groups_imp$groups_imp`),]
                                            
groups_for_pheatmap_genuses <- read.csv("heatmap_files/additional_files_hei_classification/groups_for_pheatmap_genuses.csv")
groups_for_pheatmap <- read.csv("heatmap_files/additional_files_hei_classification/groups_for_pheatmap.csv")
abunda_clr_pheatmap <- read.csv("heatmap_files/additional_files_hei_classification/abundancias_clr_para_pheatmap.csv")
rownames(abunda_clr_pheatmap) <- abunda_clr_pheatmap$Sample
abunda_clr_pheatmap$Sample <- NULL
abunda_clr_pheatmap_matrix <- as.matrix(abunda_clr_pheatmap)
colnames(abunda_clr_pheatmap_matrix)<-gsub("_"," ",colnames(abunda_clr_pheatmap_matrix))
                                            
p2 = ComplexHeatmap::HeatmapAnnotation(Importance = anno_barplot(table_imp_mikropml_rf_for_pheatmap$Overall, gp = gpar(fill = ifelse(table_imp_mikropml_rf_for_pheatmap$`groups_imp$groups_imp` == 0, "#073444", "#bae9fa"))), height = unit(3, "cm"), HEI_Group = groups_for_pheatmap_genuses$Group, col = list(HEI_Group = c("Good HEI" = "#073444", "Poor HEI" = "#bae9fa")), annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), annotation_legend_param = list(title="", labels_gp = gpar(font=2, fontsize = 12)), border = TRUE, gap = unit(1, "mm"))

p3 = ComplexHeatmap::rowAnnotation(HEI_Group = groups_for_pheatmap$Group, col = list(HEI_Group=c("Good HEI" = "#073444", "Poor HEI" = "#bae9fa")), annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), annotation_legend_param = list(title="", labels_gp = gpar(font=2, fontsize = 12)), border = TRUE, show_legend = FALSE, show_annotation_name = FALSE)

heatmap_pheatmat <- ComplexHeatmap::Heatmap(abunda_clr_pheatmap_matrix, cluster_rows = FALSE, top_annotation = p2, right_annotation = p3, col = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n = 9)), cluster_columns = FALSE, show_row_names = FALSE, heatmap_legend_param = list(title = "", legend_height = unit(4, "cm"), labels_gp = gpar(font = 2), legend_direction = "vertical"), column_names_gp = grid::gpar(fontsize = 10.5, fontface = "bold.italic"), rect_gp = gpar(col = "#073444", lwd = 1), row_split = groups_for_pheatmap$Group, row_gap = unit(2, "mm"), row_title_gp = gpar(col = c("white","white")), column_split = table_imp_mikropml_rf_for_pheatmap$`groups_imp$groups_imp`, column_title_gp = gpar(col = c("white", "white")), column_gap = unit(2, "mm"))

draw(heatmap_pheatmat,heatmap_legend_side = "left", annotation_legend_side = "bottom")
