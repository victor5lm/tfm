library(gtsummary)
library(tidyverse)
library(ggplot2)
library(gtsummary)
library(dplyr)

metadata_study_cohort <- read.csv("DATA/metadata_study_cohort.csv")

# Study cohort: population summary
trial <- datos_ias %>% select(`Age`, `Sex`, `BMI`, `Physical activity`, `Wine consumption`, `Beer consumption`, `Liquor consumption`, `Tobacco consumption`, `HEI classification`, `Highly processed food consumption`)
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
trial_val <- metadata_validation_cohort %>% select(`Age`,`Sex`,`BMI`,`HbA1c (%)`,`HbA1c IFCC (mmol/mol)`,`Glucose (mg/dl)`,`C-reactive protein (mg/dl)`,`HOMA-IR`,`Adiponectin (ug/ml)`,`HEI classification`) %>% mutate(`HOMA-IR` = factor(`HOMA-IR`, levels = c("<1.96", "1.96 to 2.99", "\u22653")))
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
        rows = label %in% "HbA1c_IFCC (mmol/mol)",
        footnote = "<u>Legend</u>:<br>\u226438: Normal values<br>
        39 to 47: Prediabetes",
        text_format = "bold") %>% 
    modify_table_styling(
        columns = label,
        rows = label %in% "Glucose (mg/dl)",
        footnote = "<u>Legend</u>:<br>\u226499: Normal values<br>
        100 to 125: Prediabetes",
        text_format = "bold") %>% 
    modify_table_styling(
        columns = label,
        rows = label %in% "C-reactive protein (mg/dl)",
        footnote = "<u>Legend</u>:<br><1: Normal values<br>
        \u22651: Abnormal values",
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
