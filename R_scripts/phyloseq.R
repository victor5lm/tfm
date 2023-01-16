[...]

processed_tb_clr <- tb_clr[!grepl("uncultured",tb_clr$Genus),]
processed_tb_clr <- processed_tb_clr[!grepl(" Incertae_Sedis",processed_tb_clr$Genus),]
tb_clr_relevant_stuff <- cbind(processed_tb_clr$SampleID,processed_tb_clr$Abundance,processed_tb_clr$Genus)
colnames(tb_clr_relevant_stuff) <- c("SampleID","Abundance","Genus")
df_spread <- as_tibble(tb_clr_relevant_stuff) %>% spread(key = Genus, value = Abundance)
