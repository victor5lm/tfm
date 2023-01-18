# Boxplot of some microbial genera whose abundances are significantly different between HPF consumption groups (Study cohort)

tabla_prueba<-tb_clr[tb_clr$Genus %in% c(" Moryella", " Muribaculaceae", " Phocea", " Prevotellaceae NK3B31 group"), ] #tb_clr is the object obtained after using psmelt on the phyloseq object.
tabla_prueba_subset<-subset(tabla_prueba,SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"))
set.seed(2019)
ggplot(data = tabla_prueba, aes(x = `HPF group`, y = Abundance, label = SampleID, group = `HPF group`, fill = `HPF group`)) +
    geom_boxplot() +
    labs(x = "Group", y = "Abundance\n") +
    facet_wrap(~ Genus, scales = "free") + 
    stat_summary(fun = "mean", geom = "point", aes(shape = "Mean"), size = 2, color = "gold", show.legend = TRUE) + geom_text_repel(aes(label=SampleID), color = ifelse(tabla_prueba$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), "white", NA), bg.color = ifelse(tabla_prueba$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), "grey30", NA), bg.r = ifelse(tabla_prueba$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), 0.15, NA), max.overlaps = Inf, point.size = 8, min.segment.length = 0, segment.color = ifelse(tabla_prueba$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), "black", NA)) +
    geom_point(color = ifelse(tabla_prueba$SampleID %in% c("PYM107","PYM202","PYM502","PYM511","PYM310","PYM402","PYM404"), "green", "grey")) +
    geom_signif(comparisons = list(c("Low HPF consumption", "High HPF consumption")), test = "wilcox.test", map_signif_level = TRUE, textsize = 5, fontface = "bold")+
    theme(plot.title = element_text(color = "#0099f8", size = 18, face = "bold", hjust = 0.5), plot.subtitle = element_text(face = "bold.italic", hjust = 0.5), strip.text = element_text(face="bold.italic", size=10), axis.title.x = element_text(size=12, face="bold", colour = "black"), axis.title.y = element_text(size=12, face="bold", colour = "black"), axis.text.x = element_text(size=10, face="bold", colour = "black"), axis.text.y = element_text(size=10, face="bold", colour = "black"), legend.title = element_text(size=11, face="bold", colour = "white"), legend.text = element_text(face="bold"), legend.position = "bottom")
        
# Boxplot of some microbial genera whose abundances are significantly different between HEI consumption groups (Study cohort)
        
        
