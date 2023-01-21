library(gtsummary)
library(tidyverse)
library(ggplot2)
library(gtsummary)
library(dplyr)

metadata_study_cohort <- read.csv("DATA/metadata_study_cohort.csv")

# Study cohort: population summary
trial <- datos_ias %>% select(`Age`, `Sex`, `BMI`, `Physical activity`, `HEI classification`, `Highly processed food consumption`)
table <- 
    tbl_summary(
        trial,
        missing = "no" # don't list missing data separately
    )%>%
    add_n() %>% # add column with total number of non-missing observations
    bold_labels() %>%
    modify_header(label ~ "**Variable**")%>% 
    modify_table_styling(
        columns = label,
        rows = label %in% c("Wine consumption","Beer consumption","Liquor consumption"),
        footnote = " <u>Legend</u>:<br>1: Never or <1 month<br>
        2: 2-4 times per week<br>
        3: Once per day<br>
        4: 1-3 time per month<br>
        5: 5-6 times per week<br>
        6: 2-3 times per day<br>
        7: Once per week",
        text_format = "bold")%>%
    modify_caption("**Summary statistics**")

table

# Validation cohort: population summary
metadata_validation_cohort <- read.csv("DATA/metadata_validation_cohort.csv")
trial_val <- metadata_validation_cohort %>% select(`Age`,`Sex`,`BMI`,`HbA1c (%)`,`HOMA-IR`,`Adiponectin (ug/ml)`,`HEI classification`) %>% mutate(`HOMA-IR` = factor(`HOMA-IR`, levels = c("<1.96", "1.96 to 2.99", "\u22653")))
table_val <- 
    tbl_summary(
        trial_val,
        missing = "no", # don't list missing data separately
        label = "Adiponectin (ug/ml)" ~ "Adiponectin (\u03BCg/ml)"
    )%>%
    add_n() %>% # add column with total number of non-missing observations
    bold_labels() %>%
    modify_header(label ~ "**Variable**") %>% 
    modify_table_styling(
        columns = label,
        rows = label %in% "HbA1c (%)",
        footnote = "<u>Legend</u>:<br><5.7: Normal values<br>
        5.7 to 6.4: Prediabetes",
        text_format = "bold") %>% 
    modify_table_styling(
        columns = label,
        rows = label %in% "HOMA-IR",
        footnote = "<u>Legend</u>:<br><1.96: Normal values<br>
        1.96 to 2.99: There might be insulin resistance<br>
        \u22653: Insulin resistance",
        text_format = "bold") %>% 
    modify_table_styling(
        columns = label,
        rows = label %in% "Adiponectin (\u03BCg/ml)",
        footnote = "<u>Legend</u>:<br><5: Abnormal values<br>
        \u22655: Normal values",
        text_format = "bold") %>%
    modify_caption("**Summary statistics**")

table_val


### DESPUÉS DE LA CORRECCIÓN DEL 20/01: NUEVAS TABLAS

# TABLA 1: COSAS COMUNES

full_data_population_tables <- read.csv("DATA/full_data_population_tables.csv") #full_data_population_tables es la tabla con toda la info, la del study cohort y la del validation cohort toda junta

colnames(full_data_population_tables) <- gsub("\\."," ",colnames(full_data_population_tables))
trial <- full_data_population_tables %>% select(`Age`, `Sex`, `BMI`, `cohort`)
table <- tbl_summary(trial,by = cohort,missing = "no")%>%
    bold_labels() %>%
    modify_header(label ~ "**Variable**")

table

# TABLA 2: HPF

full_data_population_tables <- read.csv("DATA/full_data_population_tables.csv") #full_data_population_tables es la tabla con toda la info, la del study cohort y la del validation cohort toda junta

colnames(full_data_population_tables) <- gsub("\\."," ",colnames(full_data_population_tables))
trial <- full_data_population_tables %>% select(`Age`, `Sex`, `BMI`,`Physical activity`,`Wine consumption`,`Beer consumption`,`Liquor consumption`,`Tobacco consumption`, `Highly processed food consumption`) %>% mutate(`Physical activity` = factor(`Physical activity`, levels = c("High", "Moderate", "Low")))
table <- tbl_summary(trial,by = `Highly processed food consumption`,missing = "no")%>%
    bold_labels() %>%
    modify_header(label ~ "**Variable**")%>%
    modify_spanning_header(c("stat_1", "stat_2") ~ "**Study cohort**")%>%
    modify_header(stat_2 = "**Low HPF consumption (\u2264 15 % g/day)**, N = 27")%>%
    modify_table_styling(
        columns = label,
        rows = label %in% c("Wine consumption","Beer consumption","Liquor consumption"),
        footnote = " <u>Legend</u>:<br>1: Never or <1 month<br>
        2: 2-4 times per week<br>
        3: Once per day<br>
        4: 1-3 time per month<br>
        5: 5-6 times per week<br>
        6: 2-3 times per day<br>
        7: Once per week",
        text_format = "bold")

table

# TABLA 3: HEI

full_data_population_tables <- read.csv("DATA/full_data_population_tables.csv") #full_data_population_tables es la tabla con toda la info, la del study cohort y la del validation cohort toda junta

colnames(full_data_population_tables) <- gsub("\\."," ",colnames(full_data_population_tables))
colnames(full_data_population_tables)[17] <- "HbA1c (%)"
colnames(full_data_population_tables)[19] <- "HOMA-IR"
colnames(full_data_population_tables)[21] <- "Adiponectin (ug/ml)"
full_data_population_tables$`HOMA-IR` <- gsub(">=3","\u22653",full_data_population_tables$`HOMA-IR`)
full_data_population_tables$`Adiponectin (ug/ml)` <- gsub(">=5","\u22655",full_data_population_tables$`Adiponectin (ug/ml)`)

trial <- full_data_population_tables %>% select(`Age`, `Sex`, `BMI`,`HbA1c (%)`, `HOMA-IR`, `Adiponectin (ug/ml)`, `HEI classification`, `cohort`, `HEI group`) %>% mutate(`HOMA-IR` = factor(`HOMA-IR`, levels = c("<1.96", "1.96 to 2.99", "\u22653"))) %>% mutate(`Adiponectin (ug/ml)` = factor(`Adiponectin (ug/ml)`, levels = c("<5","\u22655"))) %>% mutate(`HbA1c (%)` = factor(`HbA1c (%)`, levels = c("<5.7","5.7 to 6.4"))) %>% mutate(`HEI classification` = factor(`HEI classification`, levels = c("Excellent","Very good","Good","Acceptable","Inappropriate")))
table <- tbl_strata(data = trial, strata = cohort, .tbl_fun = ~ .x %>% tbl_summary(by = `HEI group`, missing = "no", label = "Adiponectin (ug/ml)" ~ "Adiponectin (\u03BCg/ml)"), .header = "**{strata}**, N = {n}") %>% bold_labels() %>% modify_header(stat_1_1 = "**Good HEI (\u2265 61)**, N = 34", stat_1_2 = "**Good HEI (\u2265 61)**, N = 24") %>% modify_header(label ~ "**Variable**") %>% modify_table_styling(
    columns = label,
    rows = label %in% "HbA1c (%)",
    footnote = "<u>Legend</u>:<br><5.7: Normal values<br>
        5.7 to 6.4: Prediabetes",
    text_format = "bold") %>% 
    modify_table_styling(
        columns = label,
        rows = label %in% "HOMA-IR",
        footnote = "<u>Legend</u>:<br><1.96: Normal values<br>
        1.96 to 2.99: There might be insulin resistance<br>
        \u22653: Insulin resistance",
        text_format = "bold") %>% 
    modify_table_styling(
        columns = label,
        rows = label %in% "Adiponectin (\u03BCg/ml)",
        footnote = "<u>Legend</u>:<br><5: Abnormal values<br>
        \u22655: Normal values<br><b>NA: Data not available</b></u>",
        text_format = "bold")
table
