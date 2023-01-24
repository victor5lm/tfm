# Once we have carried out the taxonomic assignment of the reads from our study cohort,
# it is time to create a phyloseq object in R that will allow us to create
# the ML classifiers based on their HPF consumption or their HEI.

# First, we import all necessary libraries
library(phyloseq)
library(tibble)
library(qiime2R)
library(tidyr)
library(readr)
library(dplyr)
library(microbiome)
library(knitr)
library(car)
library(stats)

# Now, we can create the phyloseq object:

    # ASV table
    asv <- read_qza("DATA/filtered_table.qza")

    #Taxonomía
    taxonomyq <- read_qza("DATA/taxonomy.qza")

    #Transformar tabla
    taxtableq <- taxonomyq$data %>% as.tibble() %>% separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
    taxtableq$Kingdom <- gsub("d__", "", taxtableq$Kingdom)
    taxtableq$Phylum <- gsub("p__", "", taxtableq$Phylum)
    taxtableq$Class <- gsub("c__", "", taxtableq$Class)
    taxtableq$Order <- gsub("o__", "", taxtableq$Order)
    taxtableq$Family <- gsub("f__", "", taxtableq$Family)
    taxtableq$Genus <- gsub("g__", "", taxtableq$Genus)
    taxtableq$Species <- gsub("s__", "", taxtableq$Species)

    #Árbol filogenético que generamos anteriormente
    tree <- read_qza("DATA/rooted-tree.qza")

    #Metadata
    metadata <- read_csv("DATA/metadata_study_cohort.csv")

    #Construir Objeto phyloseq
    OTUs <- otu_table(ASV$data, taxa_are_rows = TRUE)
    tree <- phy_tree(tree$data)
    TAXq <- tax_table(as.data.frame(taxtableq) %>% select_("-Confidence") %>% column_to_rownames("Feature.ID") %>% as.matrix()) 
    sample_metadata <- sample_data(metadata %>% as.data.frame())
    sample_names(sample_metadata) <- paste(metadata$SampleID)
    physeq_processed_foods <- merge_phyloseq(OTUs, tree, TAXq, sample_metadata)

# Now we have our phyloseq object, called "physeq_processed_foods"
# This object contains absolute abundances for all identified microbial taxa
# Now, we can perform the following procedures prior to diversity analyses, ML, etc:

    # First, we need to remove all taxa linked to the Eukaryota kingdom or no kingdom:
    physeq_fil <- subset_taxa(physeq_processed_foods, Kingdom != "Unassigned" & Kingdom != "Eukaryota")

    # Next, we have to agglomerate all taxonomies to the genus level, and we remove NAs:
    phy_gen <- tax_glom(physeq_fil, taxrank = "Genus", NArm = TRUE)

    # For beta-diversity calculations, relative abundances will be necessary:
    phy_gen_comp <- microbiome::transform(phy_gen, "compositional")

    # For ML classifiers, we will need the CLR-transformed relative abundances
    phy_gen_comp_clr <- microbiome::transform(phy_gen_comp, "clr")

    # For ML classifiers, the following "df_spread" variable will be highly useful
    # "df_spread" will contain all CLR-transformed abundances
    # where columns will be taxa and rows will be sample names
    tb_clr <- psmelt(phy_gen_comp_clr)
    tb_clr_relevant_stuff <- cbind(tb_clr$SampleID,tb_clr$Abundance,tb_clr$Genus)
    colnames(tb_clr_relevant_stuff) <- c("SampleID", "Abundance", "Genus")
    df_spread <- as_tibble(tb_clr_relevant_stuff) %>% spread(key = Genus, value = Abundance)

# Once we have done the previous commands,
# we can carry out the diversity analyses.
# We will use phy_gen for alpha-diversity
# and phy_gen_comp for beta-diversity.

# Please consult "R_scripts/alpha_and_beta_diversity.R" for the detailed processes.