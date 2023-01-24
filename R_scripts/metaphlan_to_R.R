# Once we got the "merged_abundance_table_gtdb.txt" file
# in "R_scripts/wgs_data.sh",
# we need to create a phyloseq object
# for the diversity analyses and ML classifier construction

library(ape)
library(phyloseq)
library(microbiome)
library(tidyr)
library(dplyr)

# First, we need this function
    metaphlanToPhyloseq <- function(
    tax,
    metadat = NULL,
    simplenames = TRUE,
    roundtointeger = FALSE,
    split="|") {
    ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
    ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
    ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
    ## if roundtointeger=TRUE, values will be rounded to the nearest integer
    xnames = rownames(tax)
    shortnames = gsub(paste0(".+\\", split), "", xnames)
    if(simplenames) {
        rownames(tax) = shortnames
    }
    if(roundtointeger) {
        tax = round(tax * 1e4)
    }
    x2 = strsplit(xnames, split=split, fixed=TRUE)
    taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow = length(x2))
    colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
    rownames(taxmat) = rownames(tax)
    for (i in 1:nrow(taxmat)){
        taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
    }
    taxmat = gsub("[a-z]__", "", taxmat)
    taxmat = phyloseq::tax_table(taxmat)
    otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
    if(is.null(metadat)){
        res = phyloseq::phyloseq(taxmat, otutab)
    }else{
        res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
    }
    return(res)
    }

# Next, we create the phyloseq object

    # We import the metadata
    metadata_val_v1 <- read.csv('DATA/metadata_validation_cohort.csv')
    # Next, we indicate the rownames
    row.names(metadata_val_v1) <- metadata_val_v1$id_voluntario
    sample.val_v1 <- sample_data(metadata_val_v1)
    sample_names(sample.val_v1) <- gsub('^', 'profiled_V1_', 
    +                                      metadata_val_v1$id_voluntario)

    # We import MetaPhlAn relative abundances table
    abundance_val_v1 <- read_csv("DATA/metaphlan_abundances_validation_cohort.csv")
    rownames(abundance_val_v1) <- abundance_val_v1$clade_name
    abundances.val_nocladename <- abundance_val_v1 %>% select(-c(clade_name))
    rownames(abundances.val_nocladename) <- abundance_val_v1$clade_name

    # Phyloseq object generation
    phyloseq_val_v1 = metaphlanToPhyloseq(abundances.val_nocladename, metadat = sample.val_v1)

    # Finally, we add to the object the phylogenetic tree
    random_tree = rtree(ntaxa(phyloseq_val_v1), rooted = T, tip.label = taxa_names(phyloseq_val_v1))
    phyloseq_val.tree <- merge_phyloseq(phyloseq_val_v1, random_tree)
    phyloseq_val.tree

# Once the phyloseq object, with the relative abundances, is created,
# we need to do the following:

    # First, we need to agglomerate taxa to the genus level and remove NAs
    # We will use this object for beta-diversity calculations
    phy_gen_val <- tax_glom(phyloseq_val.tree, taxrank = "Genus", NArm = TRUE)

    # For alpha-diversity, we will need this object (absolute abundances)
    phy_gen_val_abs <- transform_sample_counts(phy_gen_val, function(x) 1E6 * x)
                                        
    # For ML validation, we will need the CLR-transformed relative abundances
    phyloseq_val.tree_clr <- microbiome::transform(phy_gen_val, "clr")
                                        
    # For ML validation, the following "df_spread" variable will be highly useful
    # "df_spread" will contain all CLR-transformed abundances
    # where columns will be taxa and rows will be sample names
    tb_clr <- psmelt(phyloseq_val.tree_clr)
    tb_clr_relevant_stuff <- cbind(tb_clr$SampleID,tb_clr$Abundance,tb_clr$Genus)
    colnames(tb_clr_relevant_stuff) <- c("SampleID","Abundance","Genus")
    df_spread <- as_tibble(tb_clr_relevant_stuff) %>% spread(key = Genus, value = Abundance)

# Once we have done the previous commands,
# we can carry out the diversity analyses.

# Please consult "R_scripts/alpha_and_beta_diversity.R" for the detailed processes.
