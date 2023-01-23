# Script for alpha and beta diversity analyses both both cohorts grouped based on HPF consumption 
# (study cohort) and HEI (both cohorts)

# First, we need to import all necessary libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(janitor)
library(readxl)
library(phyloseq)
library(ape)
library(mia)
library(vegan)
library(rstatix)
library(scater)
library(ggsignif)
library(phyloseq)
library(tibble)
library(qiime2R)
library(readr)
library(microbiome)
library(knitr)
library(car)
library(stats)
library(cowplot)

#---ALPHA-DIVERSITY ANALYSES---

#Study cohort: HPF classification
p_alpha <- makeTreeSummarizedExperimentFromPhyloseq(phy_gen) #phy_gen is the phyloseq object with all agglomerated absolute abundances for the HPF classification
tse.list <- list("p_alpha" = p_alpha)
tse.list <- lapply(names(tse.list), function(x){
    tse.list[[x]] <- tse.list[[x]] %>%
        estimateRichness(abund_values = "counts",
                         index = "observed",
                         name = "observed") %>%
        estimateRichness(abund_values = "counts",
                         index = "chao1",
                         name = "chao1") %>%
        estimateDiversity(abund_values = "counts", 
                          index = "shannon", 
                          name = "shannon") %>%
        estimateDiversity(abund_values = "counts", 
                          index = "gini_simpson", 
                          name = "simpson")
})
names(tse.list) <- c("PYM")
df_alpha <- as.data.frame(colData(tse.list$PYM)) %>%
    select(chao1, shannon, simpson, Grupo_HPF)
counter <<- 0
graphs_alpha <-lapply(df_alpha[ ,c("chao1", "shannon", "simpson")], 
                      function(a) 
                      {counter <<- counter + 1
                      ggplot(df_alpha, aes(x = Grupo_HPF, y = a)) +
                          geom_boxplot(aes(fill = Grupo_HPF), 
                                       alpha=.5,
                                       outlier.shape = NA) +
                          geom_signif(comparisons = list(c("Low HPF consumption", "High HPF consumption")),
                                      test = "wilcox.test", map_signif_level = TRUE, textsize = 3, fontface = "bold") +
                          geom_jitter(width = 0.2,
                                      aes(colour = Grupo_HPF), size = 1.5) +
                          ylab(gsub('_', ' ', colnames(df_alpha)[counter])) + xlab(NULL) + theme_bw() + theme(axis.title.y = element_text(size=12, face="bold.italic", colour = "black"), axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=9, face="bold", colour = "black"), legend.title=element_blank(), legend.text = element_text(face = "bold"))})
ggarrange(plotlist = graphs_alpha,
          common.legend = TRUE, legend = "bottom",
          align = "hv", ncol = 3)

## Wilcoxon tests for statistical differences between both groups
wil.test_alpha_hpf <- bind_rows(pairwise_wilcox_test(df_alpha, chao1 ~ Grupo_HPF),
                      pairwise_wilcox_test(df_alpha, shannon ~ Grupo_HPF),
                      pairwise_wilcox_test(df_alpha, simpson ~ Grupo_HPF)) %>% 
adjust_pvalue(method = "BH") %>%
add_significance()

#------------
# Study cohort: HEI classification
p_hei <- makeTreeSummarizedExperimentFromPhyloseq(phy_gen) #phy_gen is the phyloseq object with all agglomerated absolute abundances for the HEI classification
tse.list_hei <- list("p_hei" = p_hei)
tse.list_hei <- lapply(names(tse.list_hei), function(x){
    tse.list_hei[[x]] <- tse.list_hei[[x]] %>%
        estimateRichness(abund_values = "counts",
                         index = "observed",
                         name = "observed") %>%
        estimateRichness(abund_values = "counts",
                         index = "chao1",
                         name = "chao1") %>%
        estimateDiversity(abund_values = "counts", 
                          index = "shannon", 
                          name = "shannon") %>%
        estimateDiversity(abund_values = "counts", 
                          index = "gini_simpson", 
                          name = "simpson")
})
names(tse.list_hei) <- c("PYM")
df_hei <- as.data.frame(colData(tse.list_hei$PYM)) %>%
    select(chao1, shannon, simpson, Grupo_HEI)
counter <<- 0
graphs_hei <-lapply(df_hei[ ,c("chao1", "shannon", "simpson")], 
                    function(a) 
                    {counter <<- counter + 1
                    ggplot(df_hei, aes(x = Grupo_HEI, y = a)) +
                        geom_boxplot(aes(fill = Grupo_HEI), 
                                     alpha=.5,
                                     outlier.shape = NA) +
                        scale_fill_brewer(palette="Paired") +
                        scale_color_brewer(palette="Paired") +
                        geom_signif(comparisons = list(c("Good HEI", "Poor HEI")),
                                    test = "wilcox.test", map_signif_level = TRUE, textsize = 3, fontface = "bold") +
                        geom_jitter(width = 0.2,
                                    aes(colour = Grupo_HEI), size = 1.5) +
                        ylab(gsub('_', ' ', colnames(df_hei)[counter])) + xlab(NULL) + theme_bw() + theme(axis.title.y = element_text(size=12, face="bold.italic", colour = "black"), axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=9, face="bold", colour = "black"), legend.title=element_blank(), legend.text = element_text(face = "bold"))})
ggarrange(plotlist = graphs_hei,
          common.legend = TRUE, legend = "bottom",
          align = "hv", ncol = 3)

## Wilcoxon tests for statistical differences between both groups
wil.test_alpha_hei <- bind_rows(pairwise_wilcox_test(df_hei, chao1 ~ Grupo_HPF),
                      pairwise_wilcox_test(df_hei, shannon ~ Grupo_HPF),
                      pairwise_wilcox_test(df_hei, simpson ~ Grupo_HPF)) %>% 
adjust_pvalue(method = "BH") %>%
add_significance()

#------------
# Validation cohort: HEI classification
p_hei_val <- makeTreeSummarizedExperimentFromPhyloseq(phy_gen_val_abs) #phy_gen_val_abs is the phyloseq object with agglomerated absolute abundances
tse.list_val <- list("p_hei_val" = p_hei_val)
tse.list_val <- lapply(names(tse.list_val), function(x){
    tse.list_val[[x]] <- tse.list_val[[x]] %>%
        estimateRichness(abund_values = "counts",
                         index = "observed",
                         name = "observed") %>%
        estimateRichness(abund_values = "counts",
                         index = "chao1",
                         name = "chao1") %>%
        estimateDiversity(abund_values = "counts", 
                          index = "shannon", 
                          name = "shannon") %>%
        estimateDiversity(abund_values = "counts", 
                          index = "gini_simpson", 
                          name = "simpson")
})
names(tse.list_val) <- c("PYM")
df_hei_val <- as.data.frame(colData(tse.list_val$PYM)) %>%
    select(chao1, shannon, simpson, group_HEI)
counter <<- 0
graphs_hei_val <-lapply(df_hei_val[ ,c("chao1", "shannon", "simpson")], 
                    function(a) 
                    {counter <<- counter + 1
                    ggplot(df_hei_val, aes(x = group_HEI, y = a)) +
                        geom_boxplot(aes(fill = group_HEI), 
                                     alpha=.5,
                                     outlier.shape = NA) +
                        scale_fill_brewer(palette="Dark2")+
                        scale_color_brewer(palette="Dark2")+
                        geom_signif(comparisons = list(c("Good HEI", "Poor HEI")),
                                    test = "wilcox.test", map_signif_level = TRUE, textsize = 3, fontface = "bold") +
                        geom_jitter(width = 0.2,
                                    aes(colour = group_HEI), size = 1.5) +
                        ylab(gsub('_', ' ', colnames(df_hei)[counter])) + xlab(NULL) + theme_bw() + theme(axis.title.y = element_text(size=12, face="bold.italic", colour = "black"), axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=9, face="bold", colour = "black"), legend.title=element_blank(), legend.text = element_text(face = "bold"))})
ggarrange(plotlist = graphs_hei_val,
          common.legend = TRUE, legend = "bottom",
          align = "hv", ncol = 3)

## Wilcoxon tests for statistical differences between both groups
wil.test_alpha_hei_val <- bind_rows(pairwise_wilcox_test(df_hei_val, chao1 ~ group_HEI),
                      pairwise_wilcox_test(df_hei_val, shannon ~ group_HEI),
                      pairwise_wilcox_test(df_hei_val, simpson ~ group_HEI)) %>% 
adjust_pvalue(method = "BH") %>%
add_significance()

#---BETA-DIVERSITY ANALYSES---

# Study cohort: HPF classification -----

## Auxiliar functions
### Calculates distances
distances <- function(study_pseq) {
    for( i in dist_methods ){
        # Calculate distance matrix
        iDist <- phyloseq::distance(study_pseq, method = i)
        dlist[[i]] <- iDist
        
        # Calculate ordination
        iPCoA  <- ordinate(study_pseq, "PCoA", distance = iDist)
        pcoa_list[[i]] <- iPCoA
    }
    to_return <- list('dlist' = dlist, 'pcoa_list' = pcoa_list)
    return(to_return)
}

# Performs PERMANOVA:
adonis_calculator <- function(dlist, study_pseq) {
    
    results.group <- lapply(names(dlist), 
                            function(x) {
                                z <- adonis2(dlist[[x]] ~ Grupo_HPF, 
                                            data = data.frame(sample_data(study_pseq)))
                                return(as.data.frame(z))
                            })
    names(results.group) <- names(dlist)
    return(results.group)
}

# Makes plots:
plotter_beta <- function(dist_methods, study_pseq, pcoa_list){
    for( i in dist_methods ){
        p <- NULL
        
        p <- plot_ordination(study_pseq, pcoa_list[[i]], 
                             color = "Grupo_HPF") + 
            geom_point(size = 2) +
            stat_ellipse(aes(group = Grupo_HPF), linetype = 2) +
            theme_bw() +
            theme(plot.title = element_text(face = 'bold', size = 16),
                  axis.title = element_text(size = 14),
                  legend.title = element_blank(),
                  axis.text = element_text(size = 12),
                  legend.text = element_text(size = 12, face = "bold"))
        
        p <- p + ggtitle(i) 
        p <- p + scale_colour_brewer(type="qual", palette="Set2")
        
        plist[[i]] <- p
    }
    
    return(plist)
}

dist_methods <- unlist(distanceMethodList)[2]
dist_methods

# to do this individually for each study:
dlist <- vector("list", length(dist_methods)) # distance matrix
names(dlist) <- dist_methods
pcoa_list <- dlist # PCoA
plist <- dlist # plots

phy_gen_comp <- microbiome::transform(phy_gen, "compositional")
phy_gen_comp.beta <- distances(phy_gen_comp)
phy_gen_comp.beta.adonis <- adonis_calculator(phy_gen_comp.beta$dlist, phy_gen_comp)
phy_gen_comp.beta.plot <- plotter_beta(dist_methods, phy_gen_comp, phy_gen_comp.beta$pcoa_list)

plist <- plotter_beta(dist_methods, phy_gen_comp, phy_gen_comp.beta$pcoa_list)
ggarrange(plotlist = plist,
          common.legend = TRUE, legend = "right")

ggarrange(plist[[1]],
          common.legend = TRUE, legend = "bottom",
          align = 'hv')

phy_gen_comp.beta.adonis # To check PERMANOVA results

# Study cohort: HEI classification -----

distances <- function(study_pseq) {
    for( i in dist_methods ){
        # Calculate distance matrix
        iDist <- distance(study_pseq, method = i)
        dlist[[i]] <- iDist
        
        # Calculate ordination
        iPCoA  <- ordinate(study_pseq, "PCoA", distance = iDist)
        pcoa_list[[i]] <- iPCoA
    }
    to_return <- list('dlist' = dlist, 'pcoa_list' = pcoa_list)
    return(to_return)
}

# Performs PERMANOVA:
adonis_calculator <- function(dlist, study_pseq) {
    
    results.group <- lapply(names(dlist), 
                            function(x) {
                                z <- adonis(dlist[[x]] ~ Grupo_HEI, 
                                            data = data.frame(sample_data(study_pseq)))
                                return(as.data.frame(z$aov.tab))
                            })
    names(results.group) <- names(dlist)
    return(results.group)
}

# Makes plots:
plotter_beta <- function(dist_methods, study_pseq, pcoa_list){
    for( i in dist_methods ){
        p <- NULL
        
        p <- plot_ordination(study_pseq, pcoa_list[[i]], 
                             color = "Grupo_HEI") +
            geom_point(size = 2) +
            stat_ellipse(aes(group = Grupo_HEI), linetype = 2) +
            theme_bw() +
            theme(plot.title = element_text(face = 'bold', size = 16),
                  axis.title = element_text(size = 14),
                  legend.title = element_blank(),
                  axis.text = element_text(size = 12),
                  legend.text = element_text(size = 12, face = "bold"))
        
        p <- p + ggtitle(i) 
        p <- p + scale_colour_brewer(type="qual", palette="Paired")
        
        plist[[i]] <- p
    }
    
    return(plist)
}

dist_methods <- unlist(distanceMethodList)[2]
dist_methods

# to do this individually for each study:
dlist <- vector("list", length(dist_methods)) # distance matrix
names(dlist) <- dist_methods
pcoa_list <- dlist # PCoA
plist <- dlist # plots

phy_gen_hei_comp <- microbiome::transform(phy_gen_hei, "compositional")
phy_gen_hei_comp.beta <- distances(phy_gen_hei_comp)
phy_gen_hei_comp.beta.adonis <- adonis_calculator(phy_gen_hei_comp.beta$dlist, phy_gen_hei_comp)
phy_gen_hei_comp.beta.plot <- plotter_beta(dist_methods, phy_gen_hei_comp, phy_gen_hei_comp.beta$pcoa_list)

plist <- plotter_beta(dist_methods, phy_gen_hei_comp, phy_gen_hei_comp.beta$pcoa_list)
ggarrange(plotlist = plist,
          common.legend = TRUE, legend = "right")

ggarrange(plist[[1]],
          common.legend = TRUE, legend = "bottom",
          align = 'hv')
                      
phy_gen_hei_comp.beta.adonis # To check PERMANOVA results
                      
# Validation cohort: HEI classification -----

distances <- function(study_pseq) {
    for(i in dist_methods){
        # Calculate distance matrix
        iDist <- phyloseq::distance(study_pseq, method = i)
        dlist[[i]] <- iDist
        
        # Calculate ordination
        iPCoA  <- ordinate(study_pseq, "PCoA", distance = iDist)
        pcoa_list[[i]] <- iPCoA
    }
    to_return <- list('dlist' = dlist, 'pcoa_list' = pcoa_list)
    return(to_return)
}

# Performs PERMANOVA:
adonis_calculator <- function(dlist, study_pseq) {
    
    results.group <- lapply(names(dlist), 
                            function(x) {
                                z <- adonis2(dlist[[x]] ~ group_HEI, 
                                             data = data.frame(sample_data(study_pseq)))
                                return(as.data.frame(z))
                            })
    names(results.group) <- names(dlist)
    return(results.group)
}

# Makes plots:
plotter_beta <- function(dist_methods, study_pseq, pcoa_list){
    for( i in dist_methods ){
        p <- NULL
        
        p <- plot_ordination(study_pseq, pcoa_list[[i]], 
                             color = "group_HEI") + 
            geom_point(size = 2) +
            stat_ellipse(aes(group = group_HEI), linetype = 2) +
            theme_bw() +
            theme(plot.title = element_text(face = 'bold', size = 16),
                  axis.title = element_text(size = 14),
                  legend.title = element_blank(),
                  axis.text = element_text(size = 12),
                  legend.text = element_text(size = 12, face = "bold"))
        
        p <- p + ggtitle(i) 
        p <- p + scale_colour_brewer(type="qual", palette="Dark2")
        
        plist[[i]] <- p
    }
    
    return(plist)
}

dist_methods <- unlist(distanceMethodList)[2]
dist_methods

# to do this individually for each study:
dlist <- vector("list", length(dist_methods)) # distance matrix
names(dlist) <- dist_methods
pcoa_list <- dlist # PCoA
plist <- dlist # plots

phy_gen_hei_comp.beta_val <- distances(phy_gen_val)
phy_gen_hei_comp.beta.adonis_val <- adonis_calculator(phy_gen_hei_comp.beta_val$dlist, phy_gen_val)
phy_gen_hei_comp.beta.plot_val <- plotter_beta(dist_methods, phy_gen_val, phy_gen_hei_comp.beta_val$pcoa_list)

plist <- plotter_beta(dist_methods, phy_gen_val, phy_gen_hei_comp.beta_val$pcoa_list)
ggarrange(plotlist = plist,
          common.legend = TRUE, legend = "right")

ggarrange(plist[[1]],
          common.legend = TRUE, legend = "bottom",
          align = 'hv')
                      
phy_gen_hei_comp.beta.adonis_val # To check PERMANOVA results
