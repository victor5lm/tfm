# Commands for the processing and taxonomic assignment of WGS data (validation cohort) -----

source ~/.bashrc

# Given all .fq.gz files, we have to check first that
# all files have been stored correctly
# by using md5sum, which calculates the 
# "digital fingerprint" of a file.

    for file in $(find -type f); do echo $file; md5sum $file > $file.md5; done
    for file in $(find -type f | grep md5); do cat $file >> md5_after_unzipping.txt ; done
    md5sum --check MD5.txt > md5sum_check_log.txt 2>&1 #Everything OK

# Several samples were sequenced more than once.
# Therefore, we need to merge those files linked to the same sample:

    # We need to unzip first our .fq.gz files.
    gunzip *.gz
    # Union of .fq files referring to the same sample.
    cat V1_1120_EKDN220029060-1A_H3NKVDSX5_L1_1.fq V1_1120_EKDN220029060-1A_H5J5KDSX5_L4_1.fq > V1_1120_1.fq
    cat V1_1120_EKDN220029060-1A_H3NKVDSX5_L1_2.fq V1_1120_EKDN220029060-1A_H5J5KDSX5_L4_2.fq > V1_1120_2.fq
    cat V1_1251_EKDN220029064-1A_H3YFWDSX5_L2_1.fq V1_1251_EKDN220029064-1A_H5JLKDSX5_L3_1.fq > V1_1251_1.fq
    cat V1_1251_EKDN220029064-1A_H3YFWDSX5_L2_2.fq V1_1251_EKDN220029064-1A_H5JLKDSX5_L3_2.fq > V1_1251_2.fq
    cat V1_1310_EKDN220029070-1A_H3NKVDSX5_L1_1.fq V1_1310_EKDN220029070-1A_H5J5KDSX5_L1_1.fq > V1_1310_1.fq
    cat V1_1310_EKDN220029070-1A_H3NKVDSX5_L1_2.fq V1_1310_EKDN220029070-1A_H5J5KDSX5_L1_2.fq > V1_1310_2.fq
    cat V1_1797_EKDN220029078-1A_H3NGLDSX5_L1_1.fq V1_1797_EKDN220029078-1A_H5J5KDSX5_L1_1.fq > V1_1797_1.fq
    cat V1_1797_EKDN220029078-1A_H3NGLDSX5_L1_2.fq V1_1797_EKDN220029078-1A_H5J5KDSX5_L1_2.fq > V1_1797_2.fq
    cat V1_1976_EKDN220029080-1A_H3NKVDSX5_L1_1.fq V1_1976_EKDN220029080-1A_H5J5KDSX5_L1_1.fq > V1_1976_1.fq
    cat V1_1976_EKDN220029080-1A_H3NKVDSX5_L1_2.fq V1_1976_EKDN220029080-1A_H5J5KDSX5_L1_2.fq > V1_1976_2.fq
    cat V1_1988_EKDN220029082-1A_H3NKVDSX5_L1_1.fq V1_1988_EKDN220029082-1A_H5J5KDSX5_L1_1.fq > V1_1988_1.fq
    cat V1_1988_EKDN220029082-1A_H3NKVDSX5_L1_2.fq V1_1988_EKDN220029082-1A_H5J5KDSX5_L1_2.fq > V1_1988_2.fq
    cat V1_1993_EKDN220029084-1A_H3NKVDSX5_L1_1.fq V1_1993_EKDN220029084-1A_H5J5KDSX5_L1_1.fq > V1_1993_1.fq
    cat V1_1993_EKDN220029084-1A_H3NKVDSX5_L1_2.fq V1_1993_EKDN220029084-1A_H5J5KDSX5_L1_2.fq > V1_1993_2.fq
    cat V1_2131_EKDN220029092-1A_H3NKVDSX5_L1_1.fq V1_2131_EKDN220029092-1A_H5J5KDSX5_L3_1.fq > V1_2131_1.fq
    cat V1_2131_EKDN220029092-1A_H3NKVDSX5_L1_2.fq V1_2131_EKDN220029092-1A_H5J5KDSX5_L3_2.fq > V1_2131_2.fq
    cat V1_2141_a_EKDN220029130-1A_H35H7DSX5_L3_1.fq V1_2141_a_EKDN220029130-1A_H5CVJDSX5_L3_1.fq > V1_2141a_1.fq
    cat V1_2141_a_EKDN220029130-1A_H35H7DSX5_L3_2.fq V1_2141_a_EKDN220029130-1A_H5CVJDSX5_L3_2.fq > V1_2141a_2.fq
    cat V1_2157_EKDN220029114-1A_H3YFWDSX5_L2_1.fq V1_2157_EKDN220029114-1A_H5JLKDSX5_L3_1.fq > V1_2157_1.fq
    cat V1_2157_EKDN220029114-1A_H3YFWDSX5_L2_2.fq V1_2157_EKDN220029114-1A_H5JLKDSX5_L3_2.fq > V1_2157_2.fq
    cat V1_2167_EKDN220029124-1A_H35H7DSX5_L3_1.fq V1_2167_EKDN220029124-1A_H5CVJDSX5_L3_1.fq > V1_2167_1.fq
    cat V1_2167_EKDN220029124-1A_H35H7DSX5_L3_2.fq V1_2167_EKDN220029124-1A_H5CVJDSX5_L3_2.fq > V1_2167_2.fq
    cat V1_2168_EKDN220029125-1A_H35H7DSX5_L4_1.fq V1_2168_EKDN220029125-1A_H5CVJDSX5_L3_1.fq > V1_2168_1.fq
    cat V1_2168_EKDN220029125-1A_H35H7DSX5_L4_2.fq V1_2168_EKDN220029125-1A_H5CVJDSX5_L3_2.fq > V1_2168_2.fq
    cat V1_2170_EKDN220029127-1A_H35H7DSX5_L3_1.fq V1_2170_EKDN220029127-1A_H3WWHDSX5_L2_1.fq > V1_2170_1.fq
    cat V1_2170_EKDN220029127-1A_H35H7DSX5_L3_2.fq V1_2170_EKDN220029127-1A_H3WWHDSX5_L2_2.fq > V1_2170_2.fq
    cat V1_2172_EKDN220029129-1A_H35H7DSX5_L3_1.fq V1_2172_EKDN220029129-1A_H5CVJDSX5_L3_1.fq > V1_2172_1.fq
    cat V1_2172_EKDN220029129-1A_H35H7DSX5_L3_2.fq V1_2172_EKDN220029129-1A_H5CVJDSX5_L3_2.fq > V1_2172_2.fq
    cat V1_614_EKDN220029036-1A_H3NKVDSX5_L1_1.fq V1_614_EKDN220029036-1A_H5J5KDSX5_L3_1.fq > V1_614_1.fq
    cat V1_614_EKDN220029036-1A_H3NKVDSX5_L1_2.fq V1_614_EKDN220029036-1A_H5J5KDSX5_L3_2.fq > V1_614_2.fq
    cat V1_650_EKDN220029038-1A_H5CW3DSX5_L2_1.fq V1_650_EKDN220029038-1A_H5J5KDSX5_L3_1.fq > V1_650_1.fq
    cat V1_650_EKDN220029038-1A_H5CW3DSX5_L2_2.fq V1_650_EKDN220029038-1A_H5J5KDSX5_L3_2.fq > V1_650_2.fq
    cat V1_668_EKDN220029040-1A_H3NKVDSX5_L1_1.fq V1_668_EKDN220029040-1A_H5J5KDSX5_L3_1.fq > V1_668_1.fq
    cat V1_668_EKDN220029040-1A_H3NKVDSX5_L1_2.fq V1_668_EKDN220029040-1A_H5J5KDSX5_L3_2.fq > V1_668_2.fq
    cat V1_752_EKDN220029041-1A_H3NKVDSX5_L1_1.fq V1_752_EKDN220029041-1A_H5J5KDSX5_L3_1.fq > V1_752_1.fq
    cat V1_752_EKDN220029041-1A_H3NKVDSX5_L1_2.fq V1_752_EKDN220029041-1A_H5J5KDSX5_L3_2.fq > V1_752_2.fq
    cat V1_760_EKDN220029042-1A_H3NKVDSX5_L1_1.fq V1_760_EKDN220029042-1A_H5J5KDSX5_L3_1.fq > V1_760_1.fq
    cat V1_760_EKDN220029042-1A_H3NKVDSX5_L1_2.fq V1_760_EKDN220029042-1A_H5J5KDSX5_L3_2.fq > V1_760_2.fq
    cat V1_783_EKDN220029044-1A_H3WWHDSX5_L4_1.fq V1_783_EKDN220029044-1A_H57WKDSX5_L3_1.fq > V1_783_1.fq
    cat V1_783_EKDN220029044-1A_H3WWHDSX5_L4_2.fq V1_783_EKDN220029044-1A_H57WKDSX5_L3_2.fq > V1_783_2.fq
    cat V1_817_EKDN220029048-1A_H3NKVDSX5_L1_1.fq V1_817_EKDN220029048-1A_H5J5KDSX5_L4_1.fq > V1_817_1.fq
    cat V1_817_EKDN220029048-1A_H3NKVDSX5_L1_2.fq V1_817_EKDN220029048-1A_H5J5KDSX5_L4_2.fq > V1_817_2.fq

# Next, we can execute metaWRAP for the 
# removal of Illumina adapters, the trimming of
# low-quality ends and the removal of host contamination.

    # Installation

        conda install -y mamba
        git clone https://github.com/bxlab/metaWRAP.git

        cd metaWRAP
        mkdir BMTAGGER_INDEX
        cd BMTAGGER_INDEX
        wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
        gunzip *fa.gz
        cat *fa > hg38.fa
        rm chr*.fa
        bmtool -d hg38.fa -o hg38.bitmask
        srprism mkindex -i hg38.fa -o hg38.srprism

        # The following line is my particular path.
        # The user should indicate here their desired path.
        export PATH=$PATH:/home/estudiantes/victor/metaWRAP/bin

        mamba create -y -n metawrap python=2.7
        conda activate metawrap

        conda config --add channels defaults
        conda config --add channels conda-forge
        conda config --add channels bioconda
        conda config --add channels ursky

        mamba install --only-deps -c ursky metawrap-mg
    
    # Execution

        # Please note this is my particular path.
        export PATH=$PATH:/home/estudiantes/victor/metaWRAP/bin

        # File names were long, so they were shortened.
        # This step can be omitted if not necessary.

        for f in *.fq; do 
            mv "$f" \
                "$(awk -F'_' '{print $1"_"$2"_"$6}' <<< $f)"
        done

        # Execution of metaWRAP

        # Please note: we need to indicate the path where the .fq files are
        # and create a new folder (here named "READ_QC")
        # where FASTQC reports will be stored for each sample.
            
        for F in RAW_READS/*_1.fq; do R=${F%_*}_2.fastq; BASE=${F##*/}; SAMPLE=${BASE%_*}; metawrap read_qc -1 $F -2 $R -t 1 -o READ_QC/$SAMPLE & done

        # Simplification of file names:

        mkdir CLEAN_READS # Final clean reads will be stored here
        for i in FINAL_READ_QC/*; do #FINAL_READ_QC is a folder where all clean reads are stored
            b=${i#*/}
            mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
            mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
        done

# Finally, for the taxonomic assignment of samples,
# we need to use MetaPhlAn4:

    # Installation

    conda create --name mpa -c conda-forge -c bioconda python=3.7 metaphlan
    conda install -c conda-forge -c bioconda metaphlan
    conda activate mpa
    metaphlan --version #MetaPhlAn version 4.0.0 (22 Sep 2022)
    
    #Installation of the bowtie2 database

    cd ~/victor # This folder should be changed to the one desired by the user.
    metaphlan --install --bowtie2db ~/victor/bowtiedb

    # Execution

    conda activate mpa
    metaphlan --version #MetaPhlAn version 4.0.0 (22 Sep 2022)

    cd ./CLEAN_READS #Or the folder where our clean reads are

    touch metaphlan_script.pl # Now, we create a perl script for the MetaPhlAn execution.
    chmod a+x metaphlan_script.pl # to allow its modification.
    vi metaphlan_script.pl # to modify the file.

    # Now we copy the following text:

        #!/usr/bin/perl

        foreach $ar (@ARGV)
        {

        @field = split (/_/, $ar);
        $name = $ar;
        $name =~ s/_1.fastq//;

        $r1 = $ar;
        $r2 = $ar;
        $r2 =~ s/1.fastq/2.fastq/;

        # pasamos ambos a metaphlan:
        system("metaphlan $r1,$r2 --bowtie2db ~/victor/bowtiedb --bowtie2out $name.bowtie2.bz2 --nproc 6 --input_type fastq -o profiled_$name.txt");
        }

    ./metaphlan_script.pl *_1.fastq > metaphlanOUT.txt 2>&1

    # Finally, we can merge all "profiled_XXX.txt" files
    # obtained after these, where each file
    # refers to each sample.

    python ~/anaconda3/envs/mpa/bin/merge_metaphlan_tables.py profiled_*.txt > merged_abundance_table.txt

# WARNING: MetaPhlan4 is based on a SGB system, which makes some genera
# not to be easily identified. For us to interpret better the results of our
# taxonomic assignment, we will transform our SGB-based profiles to
# GTDB-based profiles:

    # First, we need to execute this command where we have all "profiles".
    # The python file used here should be available after the 
    # installation of MetaPhlAn4.
    for F in $(ls); do sgb_to_gtdb_profile.py -i $F -o gtdb/$F; done

    # Now, where all new GTDB-based profiles are, we execute the following command.
    # The python file used here should be available after the 
    # installation of MetaPhlAn4.
    merge_metaphlan_tables.py profiled_* -o merged_abundance_table_gtdb.txt --gtdb_profiles

# Now we have a file called "merged_abundance_table_gtdb.txt".
# We need this file to create the phyloseq object.
# Please go to "R_scripts/metaphlan_to_R.R"

