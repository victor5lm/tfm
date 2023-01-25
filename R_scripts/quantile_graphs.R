# The following R script shows how to draw the graphs in which quantiles made
# from the validation cohort were plotted, based on several biomarkers
# related to type 2 diabetes and IR: HbA1c (%), HOMA-IR, adiponectin
# and also HEI itself. Boxplots for CLR-transformed abundances of some differentially
# abundant taxa among quantiles are plotted, where CLR-transformed abundances
# are on the y axis and quantiles are on the x axis.

library(ggplot2)
library(rstatix)
library(phyloseq)
library(dplyr)
library(ggrepel)

# QUANTILES BASED ON HOMA-IR

    # First, we need to create the quantiles
    tb_phyloseq_val_genus_quantiles_homa <- tb_clr %>% mutate(quantile = ntile(`homa-ir`,4)) #tb_clr is the table with all clr-abundances, metadata and taxonomy obtained after using psmelt
    tb_phyloseq_val_genus_quantiles_homa$quantile <- as.factor(tb_phyloseq_val_genus_quantiles_homa$quantile)
    
    # Then, via the Kruskal-Wallis test, we determine which taxa are significantly different among quantiles
    for (i in unique(tb_phyloseq_val_genus_quantiles_homa$Genus)) {out <- (kruskal_test(subset(tb_phyloseq_val_genus_quantiles_homa,Genus %in% i), Abundance ~ quantile) %>% adjust_pvalue(method = "BH") %>% add_significance()); if(out$p.adj.signif != "ns"){print(i); print(out)}}

    # After the Kruskal-Wallis test, we selected CAG-115 and Pseudoflavonifractor
    # Now, we need to select these taxa and plot the boxplot for the quantiles
    important_taxa_table_homa <- tb_phyloseq_val_genus_quantiles_homa[tb_phyloseq_val_genus_quantiles_homa$Genus %in% c("CAG-115", "Pseudoflavonifractor"), ]

    # Finally, we can plot the clr-abundances vs the quantiles of the validation cohort
    ggplot(data = important_taxa_table_homa, aes(x = quantile, y = Abundance, label = SampleID, group = quantile, fill = quantile)) +
    geom_boxplot() +
    labs(x = "Quantile", y = "Abundance\n") +
    facet_wrap(~ Genus, scales = "free") + 
    stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "orange", show.legend = TRUE) + geom_text_repel(
        color = "white",     # text color
        bg.color = "grey30", # shadow color
        bg.r = 0.15,          # shadow radius
    )+
    labs(caption = sprintf("IQR: %s\n Q1: %s\t | Q3: %s\n Min: %s\t | Max: %s", IQR(tb_phyloseq_val_genus_quantiles_homa$homa.ir), quantile(tb_phyloseq_val_genus_quantiles_homa$homa.ir, 0.25), quantile(tb_phyloseq_val_genus_quantiles_homa$homa.ir, 0.75), quantile(tb_phyloseq_val_genus_quantiles_homa$homa.ir, 0), quantile(tb_phyloseq_val_genus_quantiles_homa$homa.ir, 1)),
    ) +
    theme(
        plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), plot.caption = element_text(face = "bold.italic", size = 11, hjust = 1), strip.text = element_text(face="bold.italic", size=10), axis.title.x = element_text(size=12, face="bold", colour = "black"), axis.title.y = element_text(size=12, face="bold", colour = "black"), axis.text.x = element_text(size=10, face="bold", colour = "black"), axis.text.y = element_text(size=10, face="bold", colour = "black"), legend.title = element_text(size=11, face="bold", colour = "black"))

# QUANTILES BASED ON HbA1c (%)

    # First, we need to create the quantiles
    tb_phyloseq_val_genus_quantiles_hba1c <- tb_clr %>% mutate(quantile = ntile(`hba1c_perc`,4)) #tb_clr is the table with all clr-abundances, metadata and taxonomy obtained after using psmelt
    tb_phyloseq_val_genus_quantiles_hba1c$quantile <- as.factor(tb_phyloseq_val_genus_quantiles_hba1c$quantile)
    
    # Then, via the Kruskal-Wallis test, we determine which taxa are significantly different among quantiles
    for (i in unique(tb_phyloseq_val_genus_quantiles_hba1c$Genus)) {out <- (kruskal_test(subset(tb_phyloseq_val_genus_quantiles_hba1c,Genus %in% i), Abundance ~ quantile) %>% adjust_pvalue(method = "BH") %>% add_significance()); if(out$p.adj.signif != "ns"){print(i); print(out)}}

    # After the Kruskal-Wallis test, we selected CAG-115 and Pseudoflavonifractor
    # Now, we need to select these taxa and plot the boxplot for the quantiles
    important_taxa_table_hba1c <- tb_phyloseq_val_genus_quantiles_hba1c[tb_phyloseq_val_genus_quantiles_hba1c$Genus %in% c("Bacteroides", "Dialister"), ]

    # Finally, we can plot the clr-abundances vs the quantiles of the validation cohort
    ggplot(data = important_taxa_table_hba1c, aes(x = quantile, y = Abundance, label = SampleID, group = quantile, fill = quantile)) +
    geom_boxplot() +
    labs(x = "Quantile", y = "Abundance\n") +
    facet_wrap(~ Genus, scales = "free") + 
    stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "orange", show.legend = TRUE) + geom_text_repel(
        color = "white",     # text color
        bg.color = "grey30", # shadow color
        bg.r = 0.15,          # shadow radius
    )+
    labs(caption = sprintf("IQR: %s\n Q1: %s\t | Q3: %s\n Min: %s\t | Max: %s", IQR(tb_phyloseq_val_genus_quantiles_hba1c$hba1c_perc), quantile(tb_phyloseq_val_genus_quantiles_hba1c$hba1c_perc, 0.25), quantile(tb_phyloseq_val_genus_quantiles_hba1c$hba1c_perc, 0.75), quantile(tb_phyloseq_val_genus_quantiles_hba1c$hba1c_perc, 0), quantile(tb_phyloseq_val_genus_quantiles_hba1c$hba1c_perc, 1)),
    ) +
    theme(
        plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), plot.caption = element_text(face = "bold.italic", size = 11, hjust = 1), strip.text = element_text(face="bold.italic", size=10), axis.title.x = element_text(size=12, face="bold", colour = "black"), axis.title.y = element_text(size=12, face="bold", colour = "black"), axis.text.x = element_text(size=10, face="bold", colour = "black"), axis.text.y = element_text(size=10, face="bold", colour = "black"), legend.title = element_text(size=11, face="bold", colour = "black"))

# QUANTILES BASED ON ADIPONECTIN

    # First, we need to create the quantiles
    tb_phyloseq_val_genus_quantiles_adiponectin  <- tb_clr %>% mutate(quantile = ntile(`adiponectin_ug_ml`,4)) #tb_clr is the table with all clr-abundances, metadata and taxonomy obtained after using psmelt
    tb_phyloseq_val_genus_quantiles_adiponectin$quantile <- as.factor(tb_phyloseq_val_genus_quantiles_adiponectin$quantile)
    
    # Then, via the Kruskal-Wallis test, we determine which taxa are significantly different among quantiles
    for (i in unique(tb_phyloseq_val_genus_quantiles_adiponectin$Genus)) {out <- (kruskal_test(subset(tb_phyloseq_val_genus_quantiles_adiponectin,Genus %in% i), Abundance ~ quantile) %>% adjust_pvalue(method = "BH") %>% add_significance()); if(out$p.adj.signif != "ns"){print(i); print(out)}}

    # After the Kruskal-Wallis test, we selected CAG-115 and Pseudoflavonifractor
    # Now, we need to select these taxa and plot the boxplot for the quantiles
    important_taxa_table_adiponectin <- tb_phyloseq_val_genus_quantiles_adiponectin[tb_phyloseq_val_genus_quantiles_adiponectin$Genus %in% c("Prevotella"), ]

    # Finally, we can plot the clr-abundances vs the quantiles of the validation cohort
    ggplot(data = important_taxa_table_adiponectin, aes(x = quantile, y = Abundance, label = SampleID, group = quantile, fill = quantile)) +
    geom_boxplot() +
    labs(x = "Quantile", y = "Abundance\n") +
    facet_wrap(~ Genus, scales = "free") + 
    stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "orange", show.legend = TRUE) + geom_text_repel(
        color = "white",     # text color
        bg.color = "grey30", # shadow color
        bg.r = 0.15,          # shadow radius
    )+
    labs(caption = sprintf("IQR: %s\n Q1: %s\t | Q3: %s\n Min: %s\t | Max: %s", IQR(tb_phyloseq_val_genus_quantiles_adiponectin$adiponectin_ug_ml), quantile(tb_phyloseq_val_genus_quantiles_adiponectin$adiponectin_ug_ml, 0.25), quantile(tb_phyloseq_val_genus_quantiles_adiponectin$adiponectin_ug_ml, 0.75), quantile(tb_phyloseq_val_genus_quantiles_adiponectin$adiponectin_ug_ml, 0), quantile(tb_phyloseq_val_genus_quantiles_adiponectin$adiponectin_ug_ml, 1)),
    ) +
    theme(
        plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), plot.caption = element_text(face = "bold.italic", size = 11, hjust = 1), strip.text = element_text(face="bold.italic", size=10), axis.title.x = element_text(size=12, face="bold", colour = "black"), axis.title.y = element_text(size=12, face="bold", colour = "black"), axis.text.x = element_text(size=10, face="bold", colour = "black"), axis.text.y = element_text(size=10, face="bold", colour = "black"), legend.title = element_text(size=11, face="bold", colour = "black"))

# QUANTILES BASED ON HEI

    # First, we need to create the quantiles
    tb_phyloseq_val_genus_quantiles_hei <- tb_clr %>% mutate(quantile = ntile(`HEI`,4)) #tb_clr is the table with all clr-abundances, metadata and taxonomy obtained after using psmelt
    tb_phyloseq_val_genus_quantiles_hei$quantile <- as.factor(tb_phyloseq_val_genus_quantiles_hei$quantile)
    
    # Then, via the Kruskal-Wallis test, we determine which taxa are significantly different among quantiles
    for (i in unique(tb_phyloseq_val_genus_quantiles_hei$Genus)) {out <- (kruskal_test(subset(tb_phyloseq_val_genus_quantiles_hei,Genus %in% i), Abundance ~ quantile) %>% adjust_pvalue(method = "BH") %>% add_significance()); if(out$p.adj.signif != "ns"){print(i); print(out)}}

    # After the Kruskal-Wallis test, we selected CAG-115 and Pseudoflavonifractor
    # Now, we need to select these taxa and plot the boxplot for the quantiles
    important_taxa_table_hei <- tb_phyloseq_val_genus_quantiles_hei[tb_phyloseq_val_genus_quantiles_hei$Genus %in% c("Acetatifactor","CAG-115","Lachnospira"), ]

    # Finally, we can plot the clr-abundances vs the quantiles of the validation cohort
    ggplot(data = important_taxa_table_hei, aes(x = quantile, y = Abundance, label = SampleID, group = quantile, fill = quantile)) +
    geom_boxplot() +
    labs(x = "Quantile", y = "Abundance\n") +
    facet_wrap(~ Genus, scales = "free") + 
    stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "orange", show.legend = TRUE) + geom_text_repel(
        color = "white",     # text color
        bg.color = "grey30", # shadow color
        bg.r = 0.15,          # shadow radius
    )+
    labs(caption = sprintf("IQR: %s\n Q1: %s\t | Q3: %s\n Min: %s\t | Max: %s", IQR(tb_phyloseq_val_genus_quantiles_hei$HEI), quantile(tb_phyloseq_val_genus_quantiles_hei$HEI, 0.25), quantile(tb_phyloseq_val_genus_quantiles_hei$HEI, 0.75), quantile(tb_phyloseq_val_genus_quantiles_hei$HEI, 0), quantile(tb_phyloseq_val_genus_quantiles_hei$HEI, 1)),
    ) +
    theme(
        plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), plot.caption = element_text(face = "bold.italic", size = 11, hjust = 1), strip.text = element_text(face="bold.italic", size=10), axis.title.x = element_text(size=12, face="bold", colour = "black"), axis.title.y = element_text(size=12, face="bold", colour = "black"), axis.text.x = element_text(size=10, face="bold", colour = "black"), axis.text.y = element_text(size=10, face="bold", colour = "black"), legend.title = element_text(size=11, face="bold", colour = "black"))