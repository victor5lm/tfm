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
wil.test_alpha <- bind_rows(pairwise_wilcox_test(df_alpha, chao1 ~ Grupo_HPF),
                      pairwise_wilcox_test(df_alpha, shannon ~ Grupo_HPF),
                      pairwise_wilcox_test(df_alpha, simpson ~ Grupo_HPF)) %>% 
adjust_pvalue(method = "BH") %>%
add_significance()

# Study cohort: HEI classification

p_ias <- makeTreeSummarizedExperimentFromPhyloseq(phy_gen_ias) #phy_gen_ias is the phyloseq object with all agglomerated absolute abundances for the HEI classification
tse.list_ias <- list("p_ias" = p_ias)
tse.list_ias <- lapply(names(tse.list_ias), function(x){
    tse.list_ias[[x]] <- tse.list_ias[[x]] %>%
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
names(tse.list_ias) <- c("PYM")
df_ias <- as.data.frame(colData(tse.list_ias$PYM)) %>%
    select(chao1, shannon, simpson, Grupo_HEI)
counter <<- 0
graphs_ias <-lapply(df_ias[ ,c("chao1", "shannon", "simpson")], 
                    function(a) 
                    {counter <<- counter + 1
                    ggplot(df_ias, aes(x = Grupo_HEI, y = a)) +
                        geom_boxplot(aes(fill = Grupo_HEI), 
                                     alpha=.5,
                                     outlier.shape = NA) +
                        scale_fill_brewer(palette="Paired") +
                        scale_color_brewer(palette="Paired") +
                        geom_signif(comparisons = list(c("Good HEI", "Poor HEI")),
                                    test = "wilcox.test", map_signif_level = TRUE, textsize = 3, fontface = "bold") +
                        geom_jitter(width = 0.2,
                                    aes(colour = Grupo_HEI), size = 1.5) +
                        ylab(gsub('_', ' ', colnames(df_ias)[counter])) + xlab(NULL) + theme_bw() + theme(axis.title.y = element_text(size=12, face="bold.italic", colour = "black"), axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=9, face="bold", colour = "black"), legend.title=element_blank(), legend.text = element_text(face = "bold"))})
ggarrange(plotlist = graphs_ias,
          common.legend = TRUE, legend = "bottom",
          align = "hv", ncol = 3)

## Wilcoxon tests for statistical differences between both groups
wil.test <- bind_rows(pairwise_wilcox_test(df_ias, chao1 ~ Grupo_HPF),
                      pairwise_wilcox_test(df_ias, shannon ~ Grupo_HPF),
                      pairwise_wilcox_test(df_ias, simpson ~ Grupo_HPF)) %>% 
adjust_pvalue(method = "BH") %>%
add_significance()

# Validation cohort: HEI classification

#---BETA-DIVERSITY ANALYSES---

# Study cohort: HPF classification

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

dist_methods <- unlist(distanceMethodList)[8]
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

# Study cohort: HEI classification

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

phy_gen_ias_comp <- microbiome::transform(phy_gen_ias, "compositional")
phy_gen_ias_comp.beta <- distances(phy_gen_ias_comp)
phy_gen_ias_comp.beta.adonis <- adonis_calculator(phy_gen_ias_comp.beta$dlist, phy_gen_ias_comp)
phy_gen_ias_comp.beta.plot <- plotter_beta(dist_methods, phy_gen_ias_comp, phy_gen_ias_comp.beta$pcoa_list)

plist <- plotter_beta(dist_methods, phy_gen_ias_comp, phy_gen_ias_comp.beta$pcoa_list)
ggarrange(plotlist = plist,
          common.legend = TRUE, legend = "right")

ggarrange(plist[[1]],
          common.legend = TRUE, legend = "bottom",
          align = 'hv')

# Validation cohort: HEI classification


