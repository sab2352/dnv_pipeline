# vcf_file_size_reduction.R:  reduces vcf file sizes prior to creation of input file for denovo variant filtering pipeline.

'Usage: 
  vcf_file_size_reduction.R [--read_depth=<read_depth>] [--odds_ratio=<odds_ratio>] [--case_group=<case_group>]

  Options:
  -h --help
  --read_depth=<read_depth> setting the minimum read depth for a given variant [default: 10]
  --min_AD_alt=<min_AD_alt> minimum count of alternate alleles [default: 4]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]

' -> doc

setwd("/Users/shirazbheda/Desktop/Columbia_CPMG")

library(dplyr)
library(ggplot2)
library(data.table)
library(plyr)
library(denovolyzeR)
library("biomaRt")
require(dostats)
library(docopt)
library(here)
library(logr)

arguments <- docopt(doc, version = 'vcf_file_size_reduction.R 1.2')
read_depth <- as.numeric(arguments$read_depth)
AD_alt <- as.numeric(arguments$min_AD_alt)
pc_path <- here(paste0("dnv_filtering_pipeline/Results"))
currentDate <- format(Sys.time(), '%Y_%m_%d')

#####DUE TO FILE SIZE, WE NEED TO RUN THE TRIO ANALYSIS SCRIPT AS A LOOP#####
filenames <- list.files(here("dnv_filtering_pipeline/data"), pattern="*.csv")
filenames <- as.list(filenames)

for (i in filenames){
  trios_n <- read.csv(paste0(here("dnv_filtering_pipeline/data/"), i), stringsAsFactors = F) #XX trios, ATAV rarevar output file
  trios_n <- as.data.frame(trios_n)

#select only needed columns to reduce file sizes
  trios_n_select <- trios_n[c("VAR", "SAMPLE_1", "AD_ref_1", "AD_alt_1", "DP_1",
                                                  "GQ_1", "GT_1", "family_id", "SAMPLE_2", "AD_ref_2",
                                                  "AD_alt_2", "DP_2", "GQ_2", "GT_2", "SAMPLE_3", "AD_ref_3",
                                                  "AD_alt_3", "DP_3", "GQ_3", "GT_3", "QUAL", "MostDel", "Global_Max_AF",
                                                  "Consequence", "SYMBOL", "BIOTYPE", "HGVSc", "gnomAD_exomes_POPMAX_AF",
                                                  "gnomAD_genomes_POPMAX_AF", "hiConfDeNovo", "loConfDeNovo")]

#Filter for "protein_coding" variants only
  trios_Cavatica_protein_coding_only <- subset(trios_n_select, grepl("protein_coding", trios_n_select$BIOTYPE))
  trios_Cavatica <- trios_Cavatica_protein_coding_only

#Create numeric values of read depths
  trios_Cavatica$DP_1_1 <- as.numeric(trios_Cavatica$DP_1)
  trios_Cavatica$DP_2_1 <- as.numeric(trios_Cavatica$DP_2)
  trios_Cavatica$DP_3_1 <- as.numeric(trios_Cavatica$DP_3)

#Filter for Read Depths
  trios_Cavatica_DP10 <- subset(trios_Cavatica, as.numeric(trios_Cavatica$DP_1)>=read_depth & 
                                      as.numeric(trios_Cavatica$DP_2)>=read_depth & 
                                      as.numeric(trios_Cavatica$DP_3)>=read_depth & 
                                      !(is.na(trios_Cavatica$DP_1)) & 
                                      !(is.na(trios_Cavatica$DP_2)) & 
                                      !(is.na(trios_Cavatica$DP_3)))

#Percentage alternative allele
  trios_Cavatica_DP10$perc_AA_1 <- as.integer(trios_Cavatica_DP10$AD_alt_1)/trios_Cavatica_DP10$DP_1_1
  max(trios_Cavatica_DP10$perc_AA_1)

  trios_Cavatica_DP10$perc_AA_2 <- as.integer(trios_Cavatica_DP10$AD_alt_2)/trios_Cavatica_DP10$DP_2_1
  max(trios_Cavatica_DP10$perc_AA_2)

  trios_Cavatica_DP10$perc_AA_3 <- as.integer(trios_Cavatica_DP10$AD_alt_3)/trios_Cavatica_DP10$DP_3_1
  max(trios_Cavatica_DP10$perc_AA_3)


#tagging inherited variants
  trios_Cavatica_DP10$inherited <- ifelse((as.integer(trios_Cavatica_DP10$AD_alt_1)>=AD_alt &
                                                 trios_Cavatica_DP10$perc_AA_1>=0.2) &
                                                ((as.integer(trios_Cavatica_DP10$AD_alt_2)>=AD_alt &
                                                    trios_Cavatica_DP10$perc_AA_2>=0.2) |
                                                   (as.integer(trios_Cavatica_DP10$AD_alt_3)>=AD_alt &
                                                      trios_Cavatica_DP10$perc_AA_3>=0.2)) &
                                                trios_Cavatica_DP10$GT_1 != "'0/0'", 1, 0)

#Possible de-novo results: hiConfDeNovo
  trios_Cavatica_DP10$hiConfDeNovoY <- ifelse(trios_Cavatica_DP10$hiConfDeNovo==trios_Cavatica_DP10$SAMPLE_1, 1, 0)

#Possible de-novo results: loConfDeNovo
  trios_Cavatica_DP10$loConfDeNovoY <- ifelse(trios_Cavatica_DP10$loConfDeNovo==trios_Cavatica_DP10$SAMPLE_1, 1, 0)

#identify de-novo variants: minimum 4 alternate reads, minimum 20% alternate reads, 
#>98% reference reads for parents and high confidence de-novo based on Possible
# de-novo
  trios_Cavatica_DP10$perc_RA_2 <- trios_Cavatica_DP10$AD_ref_2/trios_Cavatica_DP10$DP_2_1
  trios_Cavatica_DP10$perc_RA_3 <- trios_Cavatica_DP10$AD_ref_3/trios_Cavatica_DP10$DP_3_1

  trios_Cavatica_DP10$de_novo <- ifelse((trios_Cavatica_DP10$AD_alt_1>=4 &
                                               trios_Cavatica_DP10$perc_AA_1>=0.2 &
                                               trios_Cavatica_DP10$perc_RA_2>=0.95 &
                                               trios_Cavatica_DP10$perc_RA_3>=0.95 &
                                               trios_Cavatica_DP10$GT_1 != "'0/0'"), 1, 0)

#Combined High confidence de-novo from Possible de-novo and manual identification
  trios_Cavatica_DP10$de_novo_inherited <- ifelse(trios_Cavatica_DP10$de_novo==1 |
                                                        trios_Cavatica_DP10$hiConfDeNovoY==1, 'de-novo',
                                                      ifelse(trios_Cavatica_DP10$inherited==1, 'inherited' ,
                                                             ifelse(trios_Cavatica_DP10$GT_1 == "'0/0'", 'ref', 0)))

#Output results as input for next step
  write.table(trios_Cavatica_DP10, sep = "\t", quote=FALSE, row.names = FALSE, 
            file=paste0(pc_path, "/", currentDate, "_dnv_filtered_tagged_minDP", read_depth, ".txt"))   
  
}
