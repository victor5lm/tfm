#ML classifiers built for the study cohort, classified based on its participants' HPF consumption.

## Necessary R packages -----
library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(viridis)
require(kernlab)

if (!require("mikropml")) {
    install.packages("mikropml")
    library(mikropml)
}

# All necessary data for classifier training -----

    ## For caret-based models

        ## classes file
        class_ml <- read.table('DATA/class_hpf_groups.tsv',
                            header = TRUE, sep = '\t', # skip = 1,
                            row.names = 1)

        ## classes for the binary classifier:
        class.bi <- factor(class_ml$class)

        ## relative abundances table:
        data_ml <- read.table('DATA/data_clr.tsv',
                            header = TRUE, sep = '\t', # skip = 1,
                            row.names = 1)

        # Preprocessing -----
        # Removal of correlated variables:
        descrCorr <- cor(data_ml)
        highCorr <- findCorrelation(descrCorr, 0.90)
        data_ml.uncor_caret <- data_ml[, -highCorr]

        #Same seed as with mikropml. 80:20 training/test set split.
        set.seed(2019)
        inTrain <- createDataPartition(class.bi,
                                    p = 4/5, list = FALSE)

        trainDescr <- data_ml.uncor_caret[inTrain, ]
        testDescr <- data_ml.uncor_caret[-inTrain, ]
        trainClass.bi <- class.bi[inTrain]
        testClass.bi <- class.bi[-inTrain]

        # sanity check:
        dim(trainDescr)[1] == length(trainClass.bi) # OK
        dim(testDescr)[1] == length(testClass.bi) # OK

    ## For mikropml-based models

        data_complete <- read.table('DATA/data_clr_hpf.tsv',
                            header = TRUE, sep = '\t', # skip = 1,
                            row.names = 1)
        # Removal of correlated variables:
        descrCorr <- cor(data_complete)
        highCorr <- findCorrelation(descrCorr, 0.90)
        data_ml.uncor_mikropml <- data_complete[, -highCorr]

## CARET-BASED CLASSIFIERS -----

## Random Forest -----

    # Hyperparameter values for grid search
    tunegrid_rf <- expand.grid(.mtry = seq(1, 70))
    ntrees <- seq(1, 100)

    set.seed(2019)
    fitControl_rf <- trainControl(
        method = "LOOCV", # LOOCV as CV method
        savePredictions = TRUE,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        verboseIter = TRUE,
        search = "grid",
    )

    originals_caret <- list() # List in which all models will be stored. The best one will be chosen from here

    for (ntree in ntrees) {
        print(ntree)
        set.seed(2019)
        fit <- train(trainDescr, trainClass.bi,
                    method = "rf",
                    tuneGrid = tunegrid_rf,
                    trControl = fitControl,
                    ntree = ntree)
        key <- toString(ntree)
        originals_caret[[key]] <- fit
    }

    # Adittional model to retrieve the best classifier out of all the ones generated in the previous loop # nolint
    get_best_model <- function(modellist) {
        acclist <- list()
        for (i in names(modellist)){
            key <- i
            acc <- modellist[[i]]$results$Accuracy
            acclist[[i]] <- max(acc)
        }
        idx_best <- names(which.max(acclist))
        return(modellist[idx_best])
    }

    # Selection of the most precise model based on its accuracy
    get_best_model(originals_caret) # ntree = 7; mtry = 4

    # Feature importances
    varImp(get_best_model(originals_caret), scale = FALSE)

## Stochastic Gradient Boosting -----

    #Hyperparameter tuning:
    gbmGrid <-  expand.grid(interaction.depth = seq(1, 30),
                            n.trees = seq(1, 100),
                            shrinkage = c(0.1, 0.9, by = 0.1),
                            n.minobsinnode = 5)

    set.seed(2019)
    fitControl_gbm <- trainControl(
        method = "LOOCV",
        savePredictions = TRUE,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        verboseIter = TRUE,
        search = "grid",
    )

    set.seed(2019)
    fit_gbm <- train(trainDescr, trainClass.bi,
                    method = "gbm",
                    tuneGrid = gbmGrid,
                    trControl = fitControl_gbm,
                    bag.fraction = 0.6)

    # Feature importances
    varImp(fit_gbm$finalModel, scale = FALSE)

## RBF-Kernel SVC -----

    #Hyperparameter tuning:
    set.seed(2019)

    # Selection of sigma = mean(sigest(taste ~ ., data = train)[-2])
    svm <- ksvm(x = as.matrix(trainDescr), y = trainClass.bi, kernel = "rbfdot", prob.model = TRUE)

    svm_grid <- data.frame(sigma = kernelf(svm)@kpar$sigma,
                        C = seq(100, 2000, by = 5))

    set.seed(2019)
    fitControl_svm <- trainControl(
        method = "LOOCV",
        savePredictions = TRUE,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        verboseIter = TRUE,
        search = "grid",
    )

    set.seed(2019)
    svm_model <- train(trainDescr, trainClass.bi,
                method = "svmRadial",
                trControl = fitControl_svm,
                tuneGrid = svm_grid)

    # Feature importances
    varImp(svm_model$finalModel, scale = FALSE)

## MIKROPML-BASED CLASSIFIERS -----

## Random Forest -----

    ntrees <- seq(1, 100)
    tuning_rf <- list(mtry = seq(1, 70))
    originals_mikropml <- list()

    for (ntree in ntrees) {
        print(ntree)
        set.seed(2019)
        results_rf <- run_ml(data_ml.uncor_mikropml, "rf",
        outcome_colname = "HPF_group",
        cross_val = caret::trainControl(method = "LOOCV"),
        training_frac = 0.80, seed = 2019,
        calculate_performance = TRUE, ntree = ntree,
        hyperparameters = tuning_rf,
        find_feature_importance = TRUE)
        key <- toString(ntree)
        originals_mikropml[[key]] <- results_rf
    }

    originals_mikropml$`17`$performance$AUC
    # 0.8833333

    # Feature importances
    varImp(originals_mikropml$`17`$trained_model, scale = FALSE)

## Extreme Gradient Boosting -----

    tuning_xgbtree <- list(nrounds = seq(1, 50, by = 5),
                        max_depth = seq(5, 20),
                        colsample_bytree = seq(0.1, 0.9, by = 0.1),
                        eta = seq(0.1, 0.5, by = 0.1),
                        gamma = 0,
                        min_child_weight = seq(1, 5),
                        subsample = seq(0.1, 0.5, by = 0.1))

    set.seed(2019)
    results_xgbtree <- run_ml(data_ml.uncor_mikropml, "xgbTree",
    outcome_colname = "HPF_group",
    cross_val = caret::trainControl(method = "LOOCV"),
    training_frac = 0.80, seed = 2019,
    calculate_performance = TRUE,
    hyperparameters = tuning_xgbtree,
    find_feature_importance = TRUE)

    # Feature importances
    varImp(results_xgbtree$trained_model, scale = FALSE)

## RBF-Kernel SVC -----

    tuning_svm <- list(sigma = kernelf(svm)@kpar$sigma, C = seq(100, 2000, by = 5))

    set.seed(2019)
    results_svm <- run_ml(data_ml.uncor_mikropml, "svmRadial",
                        outcome_colname = "HPF_group",
                        cv_times = 1,
                        kfold = 48, # 80% of the total number of samples from the study cohort
                        training_frac = 0.8, seed = 2019,
                        calculate_performance = TRUE,
                        hyperparameters = tuning_svm,
                        find_feature_importance = TRUE)
    
    # Feature importances
    varImp(results_svm$trained_model, scale = FALSE)
