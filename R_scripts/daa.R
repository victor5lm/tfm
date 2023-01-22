# This script shows all procedures for the differential abundance analysis of our cohorts
# depending on the classification of individuals

library(ggplot2)
library(DESeq2)
library(dplyr)
library(phyloseq)

# Study cohort: HPF classification -----

    ## Absolute abundances are needed for this procedure
    ## We will use the phyloseq object
    ## containing the absolute abundances

    ## First, we need to transform groups so that spaces
    ## are "_", otherwise DESeq2 will not work
    phy_gen_daa <- sample_data(phy_gen)$`Highly processed food consumption` <- gsub(" ","_",sample_data(phy_gen)$`Highly processed food consumption`)
    phygen_prev <- subset_samples(phy_gen_daa)

    ## Removal of all taxa whose total number of counts is lower or equal to 10
    phygen_prev <- prune_taxa(rowSums(otu_table(phygen_prev)) >= 10, phygen_prev)
    phygen_deseq <- phyloseq_to_deseq2(phygen_prev, ~ `Highly processed food consumption`)

    ## Execution of DESeq2 for the differential abundance analysis
    phygen_deseq <- DESeq(phygen_deseq, test = "Wald", fitType = "parametric")
    res_phygen <- results(phygen_deseq, cooksCutoff = FALSE, contrast = c("Highly processed food consumption", "Low_HPF_consumption_(<_15_%_g/day)", "High_HPF_consumption_(>_15_%_g/day)"))
    alpha <- 0.05

    ## Selection of taxa whose pvalue is below 0.05
    coldtab <- res_phygen[which(res_phygen$pvalue < alpha),]
    coldtab = cbind(as(coldtab, "data.frame"), as(tax_table(phygen_prev)[rownames(coldtab),], "matrix"))

    ## Selection of all taxa whose log2FC is between 1 and -1
    subset_coldtab <- subset(coldtab, coldtab$log2FoldChange >= 1 | coldtab$log2FoldChange <= - 1)

    ## Removal of taxa named "uncultured"
    subset_coldtab <- subset_coldtab[!grepl("uncultured", subset_coldtab$Genus), ]

    ## If needed, we can plot all interesting taxa
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    ## Genus order
    x = tapply(subset_coldtab$log2FoldChange, subset_coldtab$Genus, function(x) max(x))
    x = sort(x, TRUE)
    subset_coldtab$Genus = factor(as.character(subset_coldtab$Genus), levels=names(x))
    ggplot(subset_coldtab, aes(y = log2FoldChange, x = Genus, color = Phylum)) + geom_point(size = 4) +
        theme(axis.text.x = element_text(angle = 90, face = "bold.italic", color = "black", size = 11, hjust = 1, vjust = 0.5), axis.text.y = element_text(face = "bold", color = "black"), axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"), legend.position = "right") +
        geom_hline(aes(yintercept = -1, linetype = "dashed")) +
        geom_hline(aes(yintercept = 1, linetype = "dashed")) + guides(linetype = "none") + theme(legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold.italic", size = 11)) + scale_y_continuous(n.breaks = 10)

# Study cohort: HEI classification

    ## Absolute abundances are needed for this procedure
    ## We will use the phyloseq object
    ## containing the absolute abundances

    ## First, we need to transform groups so that spaces
    ## are "_", otherwise DESeq2 will not work
    phy_gen_daa_hei <- sample_data(phy_gen)$group_HEI <- gsub(" ", "_", sample_data(phy_gen)$group_HEI)
    phygen_prev_hei <- subset_samples(phy_gen_daa_hei)

    ## Removal of all taxa whose total number of counts is lower or equal to 10
    phygen_prev_hei <- prune_taxa(rowSums(otu_table(phygen_prev_hei)) >= 10, phygen_prev_hei)
    phygen_deseq_hei <- phyloseq_to_deseq2(phygen_prev_hei, ~ group_HEI)

    ## Execution of DESeq2 for the differential abundance analysis
    phygen_deseq_hei <- DESeq(phygen_deseq_hei, test = "Wald", fitType = "parametric")
    res_phygen_hei <- results(phygen_deseq_hei, cooksCutoff = FALSE, contrast = c("group_HEI", "Good_HEI_(>_61)", "Poor_HEI_(<_61)"))
    alpha <- 0.05
    
    ## Selection of taxa whose pvalue is below 0.05
    coldtab_phygen_hei <- res_phygen_hei[which(res_phygen_hei$pvalue < alpha), ]
    coldtab_hei = cbind(as(coldtab_phygen_hei, "data.frame"), as(tax_table(phygen_prev_hei)[rownames(coldtab_phygen_hei), ], "matrix"))

    ## Selection of all taxa whose log2FC is between 1 and -1
    subset_coldtab_hei <- subset(coldtab_hei, coldtab_hei$log2FoldChange >= 1 | coldtab_hei$log2FoldChange <= -1)
    
    ## Removal of taxa named "uncultured"
    subset_coldtab_hei <- subset_coldtab_hei[!grepl("uncultured", subset_coldtab_hei$Genus), ]

    ## If needed, we can plot all interesting taxa
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    # Genus order
    x = tapply(subset_coldtab_hei$log2FoldChange, subset_coldtab_hei$Genus, function(x) max(x))
    x = sort(x, TRUE)
    subset_coldtab_hei$Genus = factor(as.character(subset_coldtab_hei$Genus), levels = names(x))
    ggplot(subset_coldtab_hei, aes(y = log2FoldChange, x = Genus)) + geom_point(size = 4) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"), legend.position = "none", axis.text.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face="bold")) +
        geom_hline(aes(yintercept = -1, linetype = "dashed", colour = "red")) +
        geom_hline(aes(yintercept = 1, linetype = "dashed", colour = "green"))

# Validation cohort: HEI classification

    ## Absolute abundances are needed for this procedure
    ## We will use the phyloseq object
    ## containing the absolute abundances

    ## First, we need to transform groups so that spaces
    ## are "_", otherwise DESeq2 will not work
    phy_gen_daa_hei_val <- sample_data(phy_gen_val_abs)$group_HEI <- gsub(" ", "_", sample_data(phy_gen_val_abs)$group_HEI)
    phygen_prev_hei_val <- subset_samples(phy_gen_daa_hei_val)

    ## Removal of all taxa whose total number of counts is lower or equal to 10
    phygen_prev_hei_val <- prune_taxa(rowSums(otu_table(phygen_prev_hei_val)) >= 10, phygen_prev_hei_val)
    phygen_deseq_val <- phyloseq_to_deseq2(phygen_prev_hei_val, ~ group_HEI)
    
    ## Execution of DESeq2 for the differential abundance analysis
    phygen_deseq_val <- DESeq(phygen_deseq_val, test = "Wald", fitType = "parametric")
    res_phygen_val <- results(phygen_deseq_val, cooksCutoff = FALSE, contrast = c("group_HEI", "Good_HEI_(>_61)", "Poor_HEI_(<_61)"))
    alpha <- 0.05

    ## Selection of taxa whose pvalue is below 0.05
    coldtab_val <- res_phygen[which(res_phygen_val$pvalue < alpha), ]
    coldtab_val = cbind(as(coldtab_val, "data.frame"), as(tax_table(phygen_prev_hei_val)[rownames(coldtab_val), ], "matrix"))

    ## Selection of all taxa whose log2FC is between 1 and -1
    subset_coldtab_val <- subset(coldtab_val, coldtab_val$log2FoldChange >=1 | coldtab$log2FoldChange <=-1)

    ## If needed, we can plot all interesting taxa
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    # Genus order
    x = tapply(subset_coldtab_val$log2FoldChange, subset_coldtab_val$Genus, function(x) max(x))
    x = sort(x, TRUE)
    subset_coldtab_val$Genus = factor(as.character(subset_coldtab_val$Genus), levels=names(x))
    ggplot(subset_coldtab_val, aes(y=log2FoldChange, x=Genus, color = Phylum)) + geom_point(size=4) +
        theme(axis.text.x = element_text(angle=90, face = "bold.italic", color = "black", size = 11, hjust = 1, vjust = 0.5), axis.text.y = element_text(face="bold", color = "black"), axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), legend.position = "right") +
        geom_hline(aes(yintercept=-1, linetype="dashed")) +
        geom_hline(aes(yintercept=1, linetype="dashed")) + guides(linetype = "none") + theme(legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold.italic", size = 11)) + scale_y_continuous(n.breaks = 10)