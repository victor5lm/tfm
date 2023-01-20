#ML classifiers built for the study cohort, classified based on its participants' HPF consumption.

## Necessary R packages
library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(viridis)

if( !require("mikropml") ){
    install.packages("mikropml")
    library(mikropml)
}

## CARET-BASED CLASSIFIERS -----

### Random Forest -----

    ## classes file
    class.RF <- read.table('DATA/classRF_grupos.tsv',
                        header = TRUE, sep = '\t',# skip = 1,
                        row.names = 1)

    ## classes for the binary classifier:
    class.bi <- factor(class.RF$class)

    ## relative abundances table:
    data.RF <- read.table('DATA/dataRF_clr.tsv',
                        header = TRUE, sep = '\t',# skip = 1,
                        row.names = 1)

    # Preprocessing -----
    # Removal of correlated variables:
    descrCorr <- cor(data.RF)
    highCorr <- findCorrelation(descrCorr, 0.90)
    data.RF.uncor <- data.RF[ , -highCorr]

    #Same seed as with mikropml. 80:20 training/test set split.
    set.seed(2019)
    inTrain <- createDataPartition(class.bi, 
                                p = 4/5, list = FALSE)

    trainDescr <- data.RF.uncor[inTrain, ]
    testDescr <- data.RF.uncor[-inTrain, ]
    trainClass.bi <- class.bi[inTrain]
    testClass.bi <- class.bi[-inTrain]

    # sanity check:
    dim(trainDescr)[1] == length(trainClass.bi) # OK
    dim(testDescr)[1] == length(testClass.bi) # OK

    # Hyperparameter values for grid search
    tunegrid<-expand.grid(.mtry=seq(1,70))
    ntrees <- seq(1,100)

    set.seed(2019)
    fitControl <- trainControl(
        method = "LOOCV", # LOOCV as CV method
        savePredictions = TRUE,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        verboseIter = TRUE,
        search = "grid",
    )

    originals_caret <- list() # List in which all models will be stored. The best one will be chosen from here

    for (ntree in ntrees){print(ntree)
        set.seed(2019)
        fit <- train(trainDescr, trainClass.bi,
                    method = 'rf',
                    tuneGrid = tunegrid,
                    trControl = fitControl,
                    ntree = ntree)
        key <- toString(ntree)
        originals_caret[[key]] <- fit
    }

    # Adittional model to retrieve the best classifier out of all the ones generated in the previous loop
    get_best_model <- function(modellist){
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
    get_best_model(originals_caret) #mtry = 7; ntree = 4

    # 20 most important microbial taxa (or features)

### Stochastic Gradient Boosting -----

    ## classes file
    class.RF <- read.table('DATA/classRF_grupos.tsv',
                        header = TRUE, sep = '\t',# skip = 1,
                        row.names = 1)

    ## classes for the binary classifier:
    class.bi <- factor(class.RF$class)

    ## relative abundances table:
    data.RF <- read.table('DATA/dataRF_clr.tsv',
                        header = TRUE, sep = '\t',# skip = 1,
                        row.names = 1)

    # Preprocessing -----
    # Remove correlated variables:
    descrCorr <- cor(data.RF)
    highCorr <- findCorrelation(descrCorr, 0.90)
    data.RF.uncor <- data.RF[ , -highCorr]

    #mismo seed que con mikropml. Nótese que es 80:20.
    set.seed(2019)
    inTrain <- createDataPartition(class.bi, 
                                p = 4/5, list = FALSE)

    trainDescr <- data.RF.uncor[inTrain, ]
    testDescr <- data.RF.uncor[-inTrain, ]
    trainClass.bi <- class.bi[inTrain]
    testClass.bi <- class.bi[-inTrain]

    # sanity check:
    dim(trainDescr)[1] == length(trainClass.bi) #OK
    dim(testDescr)[1] == length(testClass.bi) #OK

    #Hyperparameter tuning:
    gbmGrid <-  expand.grid(interaction.depth = seq(1,30), 
                            n.trees = seq(1,100), 
                            shrinkage = c(0.1,0.9,by=0.1),
                            n.minobsinnode = 10)

    set.seed(2019)
    fitControl_gbm <- trainControl(
        method = "LOOCV",
        savePredictions = TRUE,
        classProbs = TRUE,
        summaryFunction = multiClassSummary,
        verboseIter = TRUE,
        search = "grid",
    )

    set.seed(2019)
    fit_gbm <- train(trainDescr, trainClass.bi,
                    method = 'gbm',
                    tuneGrid = gbmGrid,
                    trControl = fitControl_gbm,
                    bag.fraction = 0.7)