# In this R script, all boxplots showed in the MSc Thesis regarding significantly 
# different taxa in abundance between HPF or HEI groups are explained.
# By means of the following code, boxplots of CLR-transformed abundances of 
# significantly different taxa in abundance between groups are generated,
# where CLR-transformed abundances are plotted on the y axis and 
# the HPF/HEI groups are plotted on the x axis. One boxplot is plotted
# for each taxon, and then a final figure is obtained, which are the ones
# shown in the MSc Thesis.

library(ggplot2)
library(ggrepel)
library(ggsignif)
library(rstatix)

# Boxplot of some microbial genera whose abundances are significantly different between HPF consumption groups (Study cohort)

    # Wilcoxon rank sum test to determine those genera whose abundances are significantly different between groups
    for (i in unique(tb_clr$Genus)) {out <- (wilcox_test(subset(tb_clr,Genus %in% i), Abundance ~ HPF_group) %>% adjust_pvalue(method = "BH") %>% add_significance()); if(out$p.adj.signif != "ns"){print(i); print(out)}}

    table_taxa_hpf <- tb_clr[tb_clr$Genus %in% c(" Moryella", " Muribaculaceae", " Phocea", " Prevotellaceae NK3B31 group"), ] #tb_clr is the object obtained after using psmelt on the phyloseq object.
    
    # The following line is to highlight those individuals with T2D family history in green.
    table_taxa_hpf_subset <- subset(table_taxa_hpf, SampleID %in% c("PYM107", "PYM202", "PYM502", "PYM511", "PYM310", "PYM402", "PYM404"))
    set.seed(2019)
    ggplot(data = table_taxa_hpf, aes(x = `HPF_group`, y = Abundance, label = SampleID, group = `HPF_group`, fill = `HPF_group`)) +
        geom_boxplot() +
        labs(x = "Group", y = "Abundance\n") +
        facet_wrap(~ Genus, scales = "free") + 
        stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "gold", show.legend = TRUE) + geom_text_repel(aes(label=SampleID), color = ifelse(table_taxa_hpf$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), "white", NA), bg.color = ifelse(table_taxa_hpf$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), "grey30", NA), bg.r = ifelse(table_taxa_hpf$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), 0.15, NA), max.overlaps = Inf, point.size = 8, min.segment.length = 0, segment.color = ifelse(table_taxa_hpf$SampleID %in% c("PYM107", "PYM202", "PYM502", "PYM511", "PYM310", "PYM402", "PYM404"), "black", NA)) +
        geom_point(color = ifelse(table_taxa_hpf$SampleID %in% c("PYM107", "PYM202", "PYM502", "PYM511", "PYM310", "PYM402", "PYM404"), "green", "grey")) +
        geom_signif(comparisons = list(c("Low_HPF_consumption", "High_HPF_consumption")), test = "wilcox.test", map_signif_level = TRUE, textsize = 5, fontface = "bold")+
        theme(plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5), plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), strip.text = element_text(face = "bold.italic", size = 10), axis.title.x = element_text(size = 12, face = "bold", colour = "black"), axis.title.y = element_text(size = 12, face = "bold", colour = "black"), axis.text.x = element_text(size = 10, face = "bold", colour = "black"), axis.text.y = element_text(size = 10, face = "bold", colour = "black"), legend.title = element_text(size = 11, face = "bold", colour = "white"), legend.text = element_text(face = "bold"), legend.position = "bottom")
            
# Boxplot of some microbial genera whose abundances are significantly different between HEI consumption groups (Study cohort)

    # Wilcoxon rank sum test to determine those genera whose abundances are significantly different between groups
    for (i in unique(tb_clr$Genus)) {out <- (wilcox_test(subset(tb_clr,Genus %in% i), Abundance ~ HEI_group) %>% adjust_pvalue(method = "BH") %>% add_significance()); if(out$p.adj.signif != "ns"){print(i); print(out)}}

    table_taxa_hei <- tb_clr[tb_clr$Genus %in% c(" Butyricimonas", " CAG-352", " Candidatus_Soleaferrea", " Intestinimonas"), ] #tb_clr is the object obtained after using psmelt on the phyloseq object.
    
    # The following line is to highlight those individuals with T2D family history in green.
    table_taxa_subset_hei <- subset(table_taxa_hei,SampleID %in% c("PYM107", "PYM202", "PYM502", "PYM511", "PYM310", "PYM402", "PYM404"))
    set.seed(2019)
    ggplot(data = table_taxa_hei, aes(x = `HEI_group`, y = Abundance, label = SampleID, group = `HEI_group`, fill = `HEI_group`)) +
        geom_boxplot() +
        labs(x = "Group", y = "Abundance\n") +
        facet_wrap(~ Genus, scales = "free") + 
        stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "gold", show.legend = TRUE) + geom_text_repel(aes(label=SampleID),
        color = ifelse(table_taxa_hei$SampleID %in% c("PYM107", "PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), "white", NA),     # text color
        bg.color = ifelse(table_taxa_hei$SampleID %in% c("PYM107", "PYM202","PYM502", "PYM511", "PYM310", "PYM402", "PYM404"), "grey30", NA), # shadow color
        bg.r = ifelse(table_taxa_hei$SampleID %in% c("PYM107", "PYM202", "PYM502", "PYM511", "PYM310", "PYM402", "PYM404"), 0.15, NA),          # shadow radius
        max.overlaps = Inf,
        point.size = 8,
        min.segment.length = 0,
        segment.color = ifelse(table_taxa_hei$SampleID %in% c("PYM107", "PYM202", "PYM502", "PYM511", "PYM310", "PYM402", "PYM404"), "black", NA)) +
        geom_point(color = ifelse(table_taxa_hei$SampleID %in% c("PYM107", "PYM202", "PYM502", "PYM511", "PYM310", "PYM402", "PYM404"), "green", "grey")) +
        geom_signif(comparisons = list(c("Good_HEI", "Poor_HEI")),
                    test = "wilcox.test", map_signif_level = TRUE, textsize = 5, fontface = "bold", vjust = 0.5)+ scale_fill_brewer(palette="Dark2")+
        theme(
            plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), strip.text = element_text(face = "bold.italic", size = 10), axis.title.x = element_text(size = 12, face = "bold", colour = "black"), axis.title.y = element_text(size = 12, face = "bold", colour = "black"), axis.text.x = element_text(size = 10, face="bold", colour = "black"), axis.text.y = element_text(size = 10, face = "bold", colour = "black"), legend.title = element_text(size = 11, face = "bold", colour = "white"), legend.text = element_text(face="bold"), legend.position = "bottom")        

# Boxplot of some microbial genera whose abundances are significantly different between HEI consumption groups (Validation cohort)

    # Wilcoxon rank sum test to determine those genera whose abundances are significantly different between groups
    for (i in unique(tb_clr$Genus)) {out <- (wilcox_test(subset(tb_clr,Genus %in% i), Abundance ~ HEI_group) %>% adjust_pvalue(method = "BH") %>% add_significance()); if(out$p.adj.signif != "ns"){print(i); print(out)}}

    table_taxa_hei_val <- tb_clr[tb_clr$Genus %in% c("CAG-115", "Lachnospira"), ] #tb_clr is the object obtained after using psmelt on the phyloseq object.
    set.seed(2019)
    ggplot(data = table_taxa_hei_val, aes(x = `HEI_group`, y = Abundance, label = SampleID, group = `HEI_group`, fill = `HEI_group`)) +
        geom_boxplot() +
        labs(x = "Group", y = "Abundance\n") +
        facet_wrap(~ Genus, scales = "free") + 
        stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "gold", show.legend = TRUE) +
        geom_signif(comparisons = list(c("Good_HEI", "Poor_HEI")),
                    test = "wilcox.test", map_signif_level = TRUE, textsize = 5, fontface = "bold", vjust = 0.5)+ scale_fill_brewer(palette="Paired")+geom_point(color = "grey")+
        theme(
            plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), strip.text = element_text(face = "bold.italic", size = 10), axis.title.x = element_text(size=12, face="bold", colour = "black"), axis.title.y = element_text(size = 12, face="bold", colour = "black"), axis.text.x = element_text(size = 10, face = "bold", colour = "black"), axis.text.y = element_text(size = 10, face = "bold", colour = "black"), legend.title = element_text(size = 11, face="bold", colour = "white"), legend.text = element_text(face = "bold"), legend.position = "bottom")
