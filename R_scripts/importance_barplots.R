#This R script shows all commands for the construction of imoprtance barplots
# for all classifiers shown in the TFM document

library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(caret)
library(mikropml)
library(readr)

# HPF CLASSIFICATION -----

    # CARET-BASED MODELS -----

        # RANDOM FOREST -----

        table_imp_caret_rf <- read.csv("DATA/importances_HPF/imp_caret_rf.csv")
        imps_caret_rf <- table_imp_caret_rf %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall),
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) +
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("a") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # STOCHASTIC GRADIENT BOOSTING -----

        table_imp_caret_gbm <- read.csv("DATA/importances_HPF/imp_caret_gbm.csv")
        imps_caret_gbm <- table_imp_caret_gbm %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall),
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("b") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # RBF-KERNEL SVC -----

        table_imp_caret_svm <- read.csv("DATA/importances_HPF/imp_caret_svm.csv")
        imps_caret_svm <- table_imp_caret_svm %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("c") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

    # MIKROPML-BASED MODELS -----

        # RANDOM FOREST -----

        table_imp_mikropml_rf <- read.csv("DATA/importances_HPF/imp_mikropml_rf.csv")
        imps_mikropml_rf <- table_imp_mikropml_rf %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("d") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # EXTREME GRADIENT BOOSTING -----

        table_imp_mikropml_xgbtree <- read.csv("DATA/importances_HPF/imp_mikropml_xgbtree.csv")
        imps_mikropml_xgbtree <- table_imp_mikropml_xgbtree[1,] %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("e") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # RBF-KERNEL SVC -----

        table_imp_mikropml_svm <- read.csv("DATA/importances_HPF/imp_mikropml_svm.csv")
        imps_mikropml_svm <- table_imp_mikropml_svm %>% arrange(Overall) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("f") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

# Finally, we need to arrange all ROC plots in one picture, doing the following:

ggarrange(imps_caret_rf, imps_caret_gbm, imps_caret_svm, imps_mikropml_rf, imps_mikropml_xgbtree, imps_mikropml_svm, ncol = 2, nrow = 3, align = "hv")

# HEI CLASSIFICATION -----

    # CARET-BASED MODELS -----

        # RANDOM FOREST -----

        table_imp_caret_rf_hei <- read.csv("DATA/importances_HEI/imp_caret_rf_hei.csv")
        imps_caret_rf_hei <- table_imp_caret_rf_hei[1:10,] %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("a") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # STOCHASTIC GRADIENT BOOSTING -----

        table_imp_caret_gbm_hei <- read.csv("DATA/importances_HEI/imp_caret_gbm_hei.csv")
        imps_caret_gbm_hei <- table_imp_caret_gbm_hei[1:2,] %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("b") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # RBF-KERNEL SVC -----

        table_imp_caret_svm_hei <- read.csv("DATA/importances_HEI/imp_caret_svm_hei.csv")
        imps_caret_svm_hei <- table_imp_caret_svm_hei %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("c") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

    # MIKROPML-BASED MODELS -----

        # RANDOM FOREST -----

        table_imp_mikropml_rf_hei <- read.csv("DATA/importances_HEI/imp_mikropml_rf_hei.csv")
        imps_mikropml_rf_hei <- table_imp_mikropml_rf_hei %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("d") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # EXTREME GRADIENT BOOSTING -----

        table_imp_mikropml_xgbtree_hei <- read.csv("DATA/importances_HEI/imp_mikropml_xgbtree_hei.csv")
        imps_mikropml_xgbtree_hei <- table_imp_mikropml_xgbtree_hei[1,] %>% arrange(desc(Overall)) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = "dodge") +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = 'Features',
         fill = "Model") +
    theme_bw() + 
    ggtitle("e") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

        # RBF-KERNEL SVC -----

        table_imp_mikropml_svm_hei <- read.csv("DATA/importances_HPF/imp_caret_rf.csv")
        imps_mikropml_svm_hei <- table_imp_mikropml_svm_hei %>% arrange(Overall) %>%
    ggplot(aes(x = reorder(feature, Overall), 
               y = Overall, 
               fill = type)) +
    geom_col(width = 0.7, position = 'dodge') +
    scale_fill_manual(values = rocket(15)[7]) + # Barplot colour can be changed here, using viridis or RColorBrewer
    coord_flip() +
    labs(y = "Importance", x = "Features",
         fill = "Model") +
    theme_bw() + 
    ggtitle("f") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), axis.text.x = element_text(size = 11, face = "bold", color = "black"), legend.position = 'none') +
    theme(axis.text.y = ggtext::element_markdown())

# Finally, we need to arrange all ROC plots in one picture, doing the following:

ggarrange(imps_caret_rf_hei, imps_caret_gbm_hei, imps_caret_svm_hei, imps_mikropml_rf_hei, imps_mikropml_xgbtree_hei, imps_mikropml_svm_hei, ncol = 2, nrow = 3, align = "hv")
