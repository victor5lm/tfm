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

# Study cohort: HEI classification

# Validation cohort: HEI classification

#---BETA-DIVERSITY ANALYSES---

# Study cohort: HPF classification

# Study cohort: HEI classification

# Validation cohort: HEI classification
