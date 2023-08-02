# hard_filter_pipeline.R:  application of hard filters for denovo variant filtering pipeline.

'Usage: 
  hard_filter_pipeline.R [--read_depth=<read_depth>] [--case=<case>] 
  [--min_AD_alt=<min_AD_alt>] [--external_MAF=<external_MAF>] 
  [--binomial_test=<binomial_test>] [--prop_missing_allele=<>prop_missing_allele]
  [--geno_qual=<geno_qual>]
  
  Options:
  -h --help
  --case=<case> set the total number of samples in the analysis [default: 1]
  --read_depth=<read_depth> setting the minimum read depth for a given variant [default: 10]
  --min_AD_alt=<min_AD_alt> minimum count of alternate alleles [default: 4]
  --external_MAF=<external_MAF> threshold for external MAF filter [default: 0.00001]
  --binomial_test=<binomial_test> threshold for the binomial alternate AD test [default: 0.00001]
  --prop_missing_allele=<prop_missing_allele> threshold for max proportion of missing allele [default: .20]
  --geno_qual=<geno_qual> minimum threshold for genotype quality [default: 20]
  
' -> doc

setwd("/Users/shirazbheda/Desktop/Columbia_CPMG")

library(dplyr)
library(ggplot2)
library(data.table)
library(plyr)
library(denovolyzeR)
library("biomaRt")
require(plyr)
require(dostats)
library(docopt)
library(here)
library(logr)

arguments <- docopt(doc, version = 'hard_filter_pipeline.R 1.2')
read_depth <- as.numeric(arguments$read_depth)
min_AD_alt <- as.numeric(arguments$min_AD_alt)
external_MAF <- as.numeric(arguments$external_MAF)
binomial_test <- as.numeric(arguments$binomial_test)
prop_missing_allele <- as.numeric(arguments$prop_missing_allele)
geno_qual <- as.numeric(arguments$geno_qual)
currentDate <- format(Sys.time(), '%Y_%m_%d')

trios_dataset <- read.csv(here("dnv_filtering_pipeline/Results/allSamples_filtered_and_tagged.csv"), sep="\t",stringsAsFactors = F)

#Internal MAF - filter de-novo variants that appear in any other parent

#Create a column that counts the frequency of a variant based on varID
trios_for_internal_MAF <- trios_dataset
rm(trios_dataset)

trios_freq <- transform(trios_for_internal_MAF, freq.var = 
                              ave(seq(nrow(trios_for_internal_MAF)),
                                  VAR, FUN=length))
rm(trios_for_internal_MAF)

#split dataset into de-novo & not de-novo
trios_freq_denovo <- subset(trios_freq, (trios_freq$de_novo_inherited=="de-novo"))
trios_freq_inherited <- subset(trios_freq, (trios_freq$de_novo_inherited!="de-novo"))
rm(trios_freq)

#create a second frequency counter of de-novo variants that counts by varIDs
trios_freq_denovo_2nd_counter <- transform(trios_freq_denovo, freq.var2=ave(seq(nrow(trios_freq_denovo)), VAR, FUN=length))
trios_freq_inherited$freq.var2 <- trios_freq_inherited$freq.var
rm(trios_freq_denovo)

#recombine datasets
trios_freq_2nd_counter <- rbind(trios_freq_denovo_2nd_counter,
                                    trios_freq_inherited)
rm(trios_freq_inherited)
rm(trios_freq_denovo_2nd_counter)

#filter out any variants in which the first frequency counter does not match the second frequency counter
trios_freq_unique <- subset(trios_freq_2nd_counter,
                                (trios_freq_2nd_counter$freq.var == trios_freq_2nd_counter$freq.var2))
rm(trios_freq_2nd_counter)

table(trios_freq_unique$de_novo_inherited)
dnv_only <- subset(trios_freq_unique, 
                   trios_freq_unique$de_novo_inherited == "de-novo")
write.csv(dnv_only,
          paste0(here("dnv_filtering_pipeline/Results/"), currentDate, "_after_internal_MAF.csv"))
rm(dnv_only)

#GQ GQ_1 >= 20, DP >= 10, DP_3 >= 10 (original DP_2 & DP_3 filters were >= 20)

trios_freq_unique_DP10 <- subset(trios_freq_unique, (trios_freq_unique$GQ_1 >= geno_qual) &
                                       (trios_freq_unique$DP_2 >= read_depth) &
                                       (trios_freq_unique$DP_3 >= read_depth))
rm(trios_freq_unique)

table(trios_freq_unique_DP10$de_novo_inherited)
dnv_only <- subset(trios_freq_unique_DP10, 
                   trios_freq_unique_DP10$de_novo_inherited == "de-novo")
write.csv(dnv_only,
          paste0(here("dnv_filtering_pipeline/Results/"),currentDate,"_after_GQ_MAF.csv"))
rm(dnv_only)

#QUAL scores - SNV > 50, Indel > 300

#divide data into SNVs and Indels
trios_freq_unique_DP10_indels <- trios_freq_unique_DP10 %>% filter(grepl('ins|del|dup', HGVSc))
trios_freq_unique_DP10_SNVs <- trios_freq_unique_DP10 %>% filter(!grepl('ins|del|dup', HGVSc))
rm(trios_freq_unique_DP10)

#Do not apply any Qual score filters to Indels
trios_freq_unique_DP10_indels_QUAL300 <- trios_freq_unique_DP10_indels

#Do not apply any Qual score filters to SNVs
trios_freq_unique_DP10_SNVs_QUAL50 <- trios_freq_unique_DP10_SNVs


#Recombine subsets
trios_freq_unique_DP10_QUAL <- rbind(trios_freq_unique_DP10_indels_QUAL300, 
                                         trios_freq_unique_DP10_SNVs_QUAL50)
rm(trios_freq_unique_DP10_indels)
rm(trios_freq_unique_DP10_SNVs)
rm(trios_freq_unique_DP10_SNVs_QUAL50)
rm(trios_freq_unique_DP10_indels_QUAL300)

table(trios_freq_unique_DP10_QUAL$de_novo_inherited)
dnv_only <- subset(trios_freq_unique_DP10_QUAL, 
                   trios_freq_unique_DP10_QUAL$de_novo_inherited == "de-novo")
write.csv(dnv_only,
          paste0(here("dnv_filtering_pipeline/Results/"),currentDate, "_after_QUAL_MAF.csv"))
rm(dnv_only)


#External MAF filter - gnomAD.Exome.global_AF, gnomAD.Genome.global.AF
trios_freq_unique_DP10_QUAL_maf10e5 <- subset(trios_freq_unique_DP10_QUAL,
                                                  (Global_Max_AF<=external_MAF | is.na(Global_Max_AF)))
rm(trios_freq_unique_DP10_QUAL)

table(trios_freq_unique_DP10_QUAL_maf10e5$de_novo_inherited)
dnv_only <- subset(trios_freq_unique_DP10_QUAL_maf10e5, 
                   trios_freq_unique_DP10_QUAL_maf10e5$de_novo_inherited == "de-novo")
write.csv(dnv_only,
          paste0(here("dnv_filtering_pipeline/Results/"), currentDate, "_after_external_MAF.csv"))
rm(dnv_only)


#Indel pre-processing

#select indels only
trios_freq_unique_DP10_QUAL_maf10e5_indels <- trios_freq_unique_DP10_QUAL_maf10e5 %>%
  filter(grepl('ins|del|dup', HGVSc))

trios_freq_unique_DP10_QUAL_maf10e5_indels$vartype <- "Indel"
table(trios_freq_unique_DP10_QUAL_maf10e5_indels$de_novo_inherited)

#apply pre-processing filters
Indels_qual <- trios_freq_unique_DP10_QUAL_maf10e5_indels
rm(trios_freq_unique_DP10_QUAL_maf10e5_indels)

#add in "Percent.Alt.Read.Binomial.P" values
Indels_qual$ref_alt <- as.numeric(Indels_qual$AD_ref_1)+
  as.numeric(Indels_qual$AD_alt_1)

Indels_qual <- subset(Indels_qual, Indels_qual$ref_alt != 0)
Indels_qual <- subset(Indels_qual, Indels_qual$AD_alt_1 >= min_AD_alt)

Percent.Alt.Read.Binomial.P <- apply(Indels_qual, 1, function( x ) {
  model_binom <- binom.test( x = as.numeric( x[4] ),
                             n = as.numeric( x[48] ), 
                             p = 0.5, 
                             alternative = "two.sided", 
                             conf.level = 0.95)
  return(Percent.Alt.Read.Binomial.P = model_binom$p.value)
})

Percent.Alt.Read.Binomial.P <- as.data.frame(Percent.Alt.Read.Binomial.P)
Indels_qual <- do.call('cbind', list(Indels_qual, Percent.Alt.Read.Binomial.P))
rm(Percent.Alt.Read.Binomial.P)

Indels_qual$qualified <- ifelse((Indels_qual$Percent.Alt.Read.Binomial.P >= binomial_test),
                                "qualified", "unqualified")

Qualified_indels <- subset(Indels_qual, Indels_qual$qualified=="qualified")
rm(Indels_qual)

table(Qualified_indels$de_novo_inherited)

#Percent missing alleles
Qualified_indels$Prop_missing_allele <- (1-((Qualified_indels$AD_alt_1 + Qualified_indels$AD_ref_1)/
                                              as.numeric(Qualified_indels$DP_1)))

Qualified_indels_low_missing_alleles <- subset(Qualified_indels, (Prop_missing_allele<=prop_missing_allele))
rm(Qualified_indels)

table(Qualified_indels_low_missing_alleles$de_novo_inherited)


#SNV pre-processing

#select SNVs only
trios_freq_unique_DP10_QUAL_maf10e5_SNVs <- trios_freq_unique_DP10_QUAL_maf10e5%>%
  filter(!grepl('ins|del|dup', HGVSc))

trios_freq_unique_DP10_QUAL_maf10e5_SNVs$vartype <- "SNV"
table(trios_freq_unique_DP10_QUAL_maf10e5_SNVs$de_novo_inherited)

#apply pre-processing filters
SNVs_qual <- trios_freq_unique_DP10_QUAL_maf10e5_SNVs
rm(trios_freq_unique_DP10_QUAL_maf10e5_SNVs)


#add in "Percent.Alt.Read.Binomial.P" values
SNVs_qual$ref_alt <- as.numeric(SNVs_qual$AD_ref_1)+
  as.numeric(SNVs_qual$AD_alt_1)

SNVs_qual <- subset(SNVs_qual, SNVs_qual$ref_alt != 0)
SNVs_qual <- subset(SNVs_qual, SNVs_qual$AD_alt_1 >= min_AD_alt)

Percent.Alt.Read.Binomial.P <- apply(SNVs_qual, 1, function( x ) {
  model_binom <- binom.test( x = as.numeric( x[4] ),
                             n = as.numeric( x[48] ), 
                             p = 0.5, 
                             alternative = "two.sided", 
                             conf.level = 0.95)
  return(Percent.Alt.Read.Binomial.P = model_binom$p.value)
})

Percent.Alt.Read.Binomial.P <- as.data.frame(Percent.Alt.Read.Binomial.P)
SNVs_qual <- do.call('cbind', list(SNVs_qual, Percent.Alt.Read.Binomial.P))
rm(Percent.Alt.Read.Binomial.P)

SNVs_qual$qualified <- ifelse((SNVs_qual$Percent.Alt.Read.Binomial.P >= binomial_test),
                              "qualified", "unqualified")

Qualified_SNVs <- subset(SNVs_qual, SNVs_qual$qualified=="qualified")
rm(SNVs_qual)

table(Qualified_SNVs$de_novo_inherited)

#Percent missing alleles
Qualified_SNVs$Prop_missing_allele <- (1-((Qualified_SNVs$AD_alt_1 + Qualified_SNVs$AD_ref_1)/
                                            as.numeric(Qualified_SNVs$DP_1)))

Qualified_SNVs_low_missing_alleles <- subset(Qualified_SNVs, (Prop_missing_allele<=prop_missing_allele))
rm(Qualified_SNVs)

table(Qualified_SNVs_low_missing_alleles$de_novo_inherited)

#combine Indels & SNVs back together
Qualified_variants_low_missing_alleles <- rbind(Qualified_SNVs_low_missing_alleles, 
                                                Qualified_indels_low_missing_alleles)
rm(trios_freq_unique_DP10_QUAL_maf10e5)
rm(Qualified_SNVs_low_missing_alleles)
rm(Qualified_indels_low_missing_alleles)

table(Qualified_variants_low_missing_alleles$de_novo_inherited)

Qualified_variants_low_missing_alleles_dnv <- subset(Qualified_variants_low_missing_alleles,
                                                     Qualified_variants_low_missing_alleles$de_novo_inherited == "de-novo")

#write output to csv#
write.csv(Qualified_variants_low_missing_alleles_dnv,
          paste0(here("dnv_filtering_pipeline/Results/"), currentDate, "dnvs_post_QC.csv"))