# Commands for the processing and taxonomic assignment of 16S rRNA data (study cohort) -----

# First, we need to download the most recent version of QIIME2
# and create a conda environment -----
    wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
    conda env create -n qiime2-2022.2 --file qiime2-2022.2-py38-linux-conda.yml

    #Optional clean-up
    rm qiime2-2022.2-py38-linux-conda.yml

# From now on, we need all the .fq.gz files

# Now, we will need an additionary script to create a "manifest" file
# with information regarding where each read is located and its ID

    touch script.pl # creation of the script

    vi script.pl # vi is used to modify the script file, which is empty so far

    # now, we copy the following code:

        #!/usr/bin/perl
        use Cwd 'abs_path';
        
        open(OUT,">data.csv");

        foreach $ar (@ARGV) {

            # To extract the $name of each file:
            $name=$ar;	
            @field = split (/_/, $ar);
            $name2 = @field[0];
            
            # To get the $path:
            $path =abs_path($ar);
            
            if ($name =~ /1.fq.gz/)
            { 
                print OUT "$name2,$path,forward\n";
            }else{  
                print OUT "$name2,$path,reverse\n";
            }
        }

        close OUT;

# Now, we execute the new script to create this file, called "data.csv"
# with the aforementioned info. Please note, this file will not be
# available in this repository since the information it contains
# is not relevant.
    perl script.pl *.fq.gz

# Before calling QIIME2, we need to add the following header to data.csv:
    sed -i '1i\sample-id,absolute-filepath,direction' data.csv

# Now, we call QIIME2 and import our data:
    conda activate qiime2-2022.2

    qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path data.csv \
    --output-path demux.qza \
    --input-format PairedEndFastqManifestPhred33 \

    qiime demux summarize \
    --i-data demux.qza \
    --o-visualization demux.qzv

# Next, we carry out the quality filtering procedure with DADA2. 
# No additional reads will be removed since the original raw reads
# already had a good quality

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table table.qza \ # Info as for number of features, number of samples, frequency per sample, etc
  --o-representative-sequences rep-seqs.qza \ # Identified features and their sequence
  --o-denoising-stats denoising-stats.qza # Statistics regarding the filtering processes

# We can then transform these .qza artefacts to visualisation files (.qzv) in case we are interested in viewing their content
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
  
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# We then perform the taxonomic assignment of our features
wget https://data.qiime2.org/2022.2/common/silva-138-99-nb-classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
  
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Next, we need to remove features linked to mitochondria or chloroplasts from table.qza and rep-seqs.qza
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table filtered_table.qza
  
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences filtered_rep-seqs.qza

# Finally, we perform the phylogenetic tree construction with MAFFT and FastTree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences filtered_rep-seqs.qza \
  --o-alignment aligned_filtered-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# We are done with the taxonomic assignment.
# Now, we will need taxonomy.qza, filtered_table.qza and rooted-tree.qza
# for the phyloseq object construction, along with the metadata.
# Please go to R_scripts/qiime_to_R.R
