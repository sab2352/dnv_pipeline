#\!/bin/bash
# 2023-03-23
# Run denovo variant pipeline 


################ Must be defined by user ###################
PROJECT="." # path to your project
# source /nfs/goldstein/software/centos7/python-3.9.7-x86_64_shared/python3.9.7-ENV.sh

case_path=$(cat ./Input/dnv_pipeline_inputs.yaml | shyaml get-value USER_VARIABLE.case_path)
  echo "case_path path is "
  echo $case_path
  echo "The current working directory is"
  pwd
  echo "These should be the same" #in the future, should error if not the same
  sleep 2

################ Other paths ###################
#atav="/nfs/goldstein/software/sh/atav.sh --email" # remove --email if you don't want to get them
#source /usr/local/igm/non-atav-tools/R-4.1.0_with_gcc_10-x86_64/R-4.1.0_with_gcc_10-x86_64-ENV.sh
# Rscript_var="/usr/local/igm/non-atav-tools/R-4.1.0_with_gcc_10-x86_64/bin/Rscript"
Rscript_var="Rscript"

################ Path to variables to be changed by user ###################
sample_vcf_path="./data/simulated_trio_data.csv"
case_group=$(cat ./Input/dnv_pipeline_inputs.yaml | shyaml get-value USER_VARIABLE.case_group)
  echo "case_group path is "
  echo $case_group

################ Paths to R scripts ###################
vcf_file_size_reduction_var=$PROJECT/Scripts/vcf_file_size_reduction.R
denovo_hard_filters_var=$PROJECT/Scripts/hard_filter_pipeline.R

if [[ -z $1 ]]
then
  echo "needs a parameter to determine which step to run"
  exit
fi

################ Read Me ###################
# Each step is run separately and then commented out so that the next step can be run
# All R code is run on dev1 to dev4
# All ATAV commands are run on qs1


########################################################################################################
########################################################################################################
################ run vcf_file_size_reduction.R to reduce file sizes and tag variants ###################
########################################################################################################
########################################################################################################
read_depth=$(cat ./Input/dnv_pipeline_inputs.yaml | shyaml get-value DNV_FILTERING_VARIABLE.read_depth)

if [[ $1 = "vcfReduction" ]]
then
  echo "Reducing vcf file sizes..."

echo "Read depth must be greater than or equal to $read_depth"
sleep 2

  # Learn about the input options for vcf_file_size_reduction.R
  $Rscript_var $vcf_file_size_reduction_var --help

  # run vcf reduction script
  $Rscript_var $vcf_file_size_reduction_var --read_depth $read_depth # --sample_path_w_filename $sample_vcf_path --date $date

fi

################ File concatenation ###################
currentDate=$(date '+%Y_%m_%d')

echo "Today's date is..."
echo $currentDate

if [[ $1 == "concatenate" ]]
then
  echo "concatenating reduced vcf files..."
  head -n1 ./Results/${currentDate}_dnv_filtered_tagged_minDP${read_depth}.txt > ./Results/allSamples_filtered_and_tagged.csv #specify a single input & output file name
  for fname in ./Results/${currentDate}_dnv_filtered_tagged_minDP${read_depth}.txt; do     tail -n+2 $fname >> ./Results/allSamples_filtered_and_tagged.csv; done
fi

########################################################################################################
########################################################################################################
###############                run hard_filter_pipeline.R                  #############################
########################################################################################################
########################################################################################################

if [[ $1 = "dnvFiltering" ]]
then
	echo "Applying filters to denovo and inherited variants..."

# read_depth=10 #this is defined by user
read_depth=$(cat ./Input/dnv_pipeline_inputs.yaml | shyaml get-value DNV_FILTERING_VARIABLE.read_depth)
#sleep 2

  # Learn about the input options for hard_filter_pipeline.R
  $Rscript_var $denovo_hard_filters_var --help
  
  # run dnv filtering pipeline
  $Rscript_var $denovo_hard_filters_var
fi


