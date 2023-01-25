# This R script shows all commands for the construction of ROC curves
# for all classifiers shown in the TFM document.

# Models used here for the ROC curve generation
# are the ones detailed in R_scripts/ML_models_HEI.R
# and R_scripts/ML_models_HPF.R.

library(ggplot2)
library(pROC)
library(dplyr)
library(caret)
library(mikropml)
library(viridis)
library(RColorBrewer)
library(purrr)
library(ggpubr)

# Helper function for ROC curve plotting:
    roc.plot <- function(roc.list, title.in, p.ci = FALSE,
                        shuffle = TRUE) {
    # generates palette
    if (shuffle == TRUE){
        lengthroc <- length(roc.list)/2
        pal <- rep(viridis::viridis(lengthroc), 2)
        aes_ggroc <- c("linetype", "colour")
        line_types <- c(rep("solid", lengthroc), rep("twodash", lengthroc))
    }
    else {
        pal <- viridis::viridis(length(roc.list)) # Plot colour is changed here
        aes_ggroc <- "colour"
        line_types <- "solid"
    }
    
    # extracts auc
    roc.list %>% 
        map(~tibble(AUC = .x$auc)) %>% 
        bind_rows(.id = "name") -> data.auc
    
    # generates labels
    data.auc %>% 
        mutate(label_long = paste0(name,", AUC = ", paste(round(AUC, 2))),
            label_AUC = paste0("AUC = ",paste(round(AUC, 2)))) -> data.labels
    
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

# HPF CLASSIFICATION -----

    # CARET-BASED MODELS -----

        # RANDOM FOREST -----

        pred_caret_rf <- predict(get_best_model(originals_caret), testDescr, type = "prob")
        roc_caret_rf <- roc(testClass.bi, pred_caret_rf$"Low_HPF_consumption")
        my.roc_caret_rf <- list("RF (caret)" = roc_caret_rf)
        plot_roc_caret_rf <- roc.plot(my.roc_caret_rf, "", p.ci = TRUE, shuffle = FALSE)

        # STOCHASTIC GRADIENT BOOSTING -----

        pred_caret_gbm <- predict(fit_gbm, testDescr, type = "prob")
        roc_caret_gbm <- roc(testClass.bi, pred_caret_gbm$"Low_HPF_consumption")
        my.roc_caret_gbm <- list("GBM (caret)" = roc_caret_gbm)
        plot_roc_caret_gbm <- roc.plot(my.roc_caret_gbm, "", p.ci = TRUE, shuffle = FALSE)

        # RBF-KERNEL SVC -----

        pred_caret_svm <- predict(svm_model, testDescr, type = "prob")
        roc_caret_svm <- roc(testClass.bi, pred_caret_svm$"Low_HPF_consumption")
        my.roc_caret_svm <- list("SVC (caret)" = roc_caret_svm)
        plot_roc_caret_svm <- roc.plot(my.roc_caret_svm, "", p.ci = TRUE, shuffle = FALSE)

    # MIKROPML-BASED MODELS -----

        # RANDOM FOREST -----

        pred_mikropml_rf <- predict(originals_mikropml$`17`$trained_model, originals_mikropml$`17`$test_data, type = "prob")
        roc_mikropml_rf <- roc(as.factor(originals_mikropml$`17`$test_data$HPF_group), pred_mikropml_rf$"Low_HPF_consumption")
        my.roc_mikropml_rf <- list("RF (mikropml)" = roc_mikropml_rf)
        plot_roc_mikropml_rf <- roc.plot(my.roc_mikropml_rf, "", p.ci = TRUE, shuffle = FALSE)

        # EXTREME GRADIENT BOOSTING -----

        pred_mikropml_xgbtree <- predict(results_xgbtree$trained_model, results_xgbtree$test_data, type = "prob")
        roc_mikropml_xgbtree <- roc(as.factor(results_xgbtree$test_data$HPF_group), pred_mikropml_xgbtree$"Low_HPF_consumption")
        my.roc_mikropml_xgbtree <- list("XGBTREE (mikropml)" = roc_mikropml_xgbtree)
        plot_roc_mikropml_xgbtree <- roc.plot(my.roc_mikropml_xgbtree, "", p.ci = TRUE, shuffle = FALSE)

        # RBF-KERNEL SVC -----

        pred_mikropml_svm <- predict(results_svm$trained_model, results_svm$test_data, type = "prob")
        roc_mikropml_svm <- roc(as.factor(results_svm$test_data$HPF_group), pred_mikropml_svm$"Low_HPF_consumption")
        my.roc_mikropml_svm <- list("SVC (mikropml)" = roc_mikropml_svm)
        plot_roc_mikropml_svm <- roc.plot(my.roc_mikropml_svm, "", p.ci = TRUE, shuffle = FALSE)

# Finally, we need to arrange all ROC plots in one picture, doing the following:

ggarrange(plot_roc_caret_rf, plot_roc_mikropml_rf, plot_roc_caret_gbm, plot_roc_mikropml_xgbtree, plot_roc_caret_svm, plot_roc_mikropml_svm)

# HEI CLASSIFICATION -----

    # CARET-BASED MODELS -----
    
        # RANDOM FOREST -----

        pred_caret_rf_hei <- predict(get_best_model(originals_caret_hei), testDescr, type = "prob")
        roc_caret_rf_hei <- roc(testClass.bi, pred_caret_rf_hei$"Good_HEI")
        my.roc_caret_rf_hei <- list("RF (caret)" = roc_caret_rf_hei)
        plot_roc_caret_rf_hei <- roc.plot(my.roc_caret_rf_hei, "", p.ci = TRUE, shuffle = FALSE)

        # STOCHASTIC GRADIENT BOOSTING -----

        pred_caret_gbm_hei <- predict(fit_gbm_hei, testDescr, type = "prob")
        roc_caret_gbm_hei <- roc(testClass.bi, pred_caret_gbm_hei$"Good_HEI")
        my.roc_caret_gbm_hei <- list("GBM (caret)" = roc_caret_gbm_hei)
        plot_roc_caret_gbm_hei <- roc.plot(my.roc_caret_gbm_hei, "", p.ci = TRUE, shuffle = FALSE)

        # RBF-KERNEL SVC -----

        pred_caret_svm_hei <- predict(svm_model_hei, testDescr, type = "prob")
        roc_caret_svm_hei <- roc(testClass.bi, pred_caret_svm_hei$"Good_HEI")
        my.roc_caret_svm_hei <- list("SVC (caret)" = roc_caret_svm_hei)
        plot_roc_caret_svm_hei <- roc.plot(my.roc_caret_svm_hei, "", p.ci = TRUE, shuffle = FALSE)

    # MIKROPML-BASED MODELS -----

        # RANDOM FOREST -----

        pred_mikropml_rf_hei <- predict(originals_mikropml$`88`$trained_model, originals_mikropml$`88`$test_data, type = "prob")
        roc_mikropml_rf_hei <- roc(as.factor(originals_mikropml$`88`$test_data$HEI_group), pred_mikropml_rf_hei$"Good_HEI")
        my.roc_mikropml_rf_hei <- list("RF (mikropml)" = roc_mikropml_rf_hei)
        plot_roc_mikropml_rf_hei <- roc.plot(my.roc_mikropml_rf_hei, "", p.ci = TRUE, shuffle = FALSE)

        # EXTREME GRADIENT BOOSTING -----

        pred_mikropml_xgbtree_hei <- predict(results_xgbtree_hei$trained_model, results_xgbtree_hei$test_data, type = "prob")
        roc_mikropml_xgbtree_hei <- roc(as.factor(results_xgbtree_hei$test_data$HEI_group), pred_mikropml_xgbtree_hei$"Good_HEI")
        my.roc_mikropml_xgbtree_hei <- list("XGBTREE (mikropml)" = roc_mikropml_xgbtree_hei)
        plot_roc_mikropml_xgbtree_hei <- roc.plot(my.roc_mikropml_xgbtree_hei, "", p.ci = TRUE, shuffle = FALSE)

        # RBF-KERNEL SVC -----

        pred_mikropml_svm_hei <- predict(results_svm_hei$trained_model, results_svm_hei$test_data, type = "prob")
        roc_mikropml_svm_hei <- roc(as.factor(results_svm_hei$test_data$HEI_group), pred_mikropml_svm_hei$"Good_HEI")
        my.roc_mikropml_svm_hei <- list("SVC (mikropml)" = roc_mikropml_svm_hei)
        plot_roc_mikropml_svm_hei <- roc.plot(my.roc_mikropml_svm_hei, "", p.ci = TRUE, shuffle = FALSE)

# Finally, we need to arrange all ROC plots in one picture, doing the following:

ggarrange(plot_roc_caret_rf_hei, plot_roc_mikropml_rf_hei, plot_roc_caret_gbm_hei, plot_roc_mikropml_xgbtree_hei, plot_roc_caret_svm_hei, plot_roc_mikropml_svm_hei)