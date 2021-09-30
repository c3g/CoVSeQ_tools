#!/bin/bash
set -eu -o pipefail

# LSPQ NANOPORE ANALYSIS
##################################################
# NANOPORE DATA TRANSFER
# Transfers output from the GenPipes Nanopore covseq pipeline into the format expected by LSPQ
##################################################

####################
if [ "$#" != 2 ]
then
  echo "Error: you provided an incorrect number of arguments." 
  echo "Usage: $(basename $0) <input_dir> <output_dir>"
  exit 1
fi

####################
# Set file paths 
INPUT_NANOPORE_PATH=${1}
OUTPUT_REPORT_PATH=${2}
RUN_NAME=$(basename ${INPUT_NANOPORE_PATH})

####################
# Create output directory
TOP_OUTPUT=${OUTPUT_REPORT_PATH}/${RUN_NAME}
ANALYSIS_OUTPUT=${OUTPUT_REPORT_PATH}/${RUN_NAME}/analysis/${RUN_NAME}_nanopolish_800x
mkdir -p ${ANALYSIS_OUTPUT}

####################
# Sync raw data
RAW_DATA_DIR=$(find $INPUT_NANOPORE_PATH -type d -name "$(date +%Y)"*)
rsync -rvlt ${RAW_DATA_DIR} ${TOP_OUTPUT}/

####################
# If basecall and demultiplex exist, sync them 
if [ $(find ${INPUT_NANOPORE_PATH} -name basecall -type d) ] ; then
    rsync -rvlt $(find ${INPUT_NANOPORE_PATH} -name basecall -type d) ${TOP_OUTPUT}
fi 

if [ $(find ${INPUT_NANOPORE_PATH} -name demultiplex -type d) ] ; then
    rsync -rvlt $(find ${INPUT_NANOPORE_PATH} -name demultiplex -type d) ${TOP_OUTPUT}/basecall/
fi     

####################
# If Genpipes analyses ouptut directories exist, sync them
for gp_out_dir in {alignment,artic_nanopolish,consensus,host_removal,variant,dna}; do 
    if [ $(find ${INPUT_NANOPORE_PATH} -name $gp_out_dir -type d) ] ; then
        rsync -rvLt $(find ${INPUT_NANOPORE_PATH} -name $gp_out_dir -type d)/* ${ANALYSIS_OUTPUT}/
    fi
done

####################
# If ncov_tools report exists, sync it 
if [ $(find ${INPUT_NANOPORE_PATH} -name ncov_tools -type d) ] ; then
    rsync -rvlt $(find ${INPUT_NANOPORE_PATH} -name ncov_tools -type d) ${ANALYSIS_OUTPUT}/
fi     

####################
# Create links for files to ensure background compatibility 
if [ $(find ${INPUT_NANOPORE_PATH} -name "readset.txt" ) ] ; then  
    # Save Readset file location
    READSET=$(find ${INPUT_NANOPORE_PATH} -name "readset.txt" )
    # Sync readset file
    cp ${READSET} ${ANALYSIS_OUTPUT}/
    # Loop through samples in readset
    for sample in $(grep -v "^Sample" ${READSET} | cut -f 1 ); do 
       echo "Linking files for ${sample}..." 
       # Find and link metrics file
       if [ $(find ${ANALYSIS_OUTPUT}/${sample}/ -name "${sample}.metrics.csv" -type f) ] ; then 
            ln -sf general_metrics/${sample}.metrics.csv ${ANALYSIS_OUTPUT}/${sample}/${sample}.metrics.csv
       fi 
       # Find and link clean fastq file
       if [  $(find ${ANALYSIS_OUTPUT}/${sample}/ -name "${sample}.host_removed.fastq.gz" -type f) ] ; then 
            ln -sf ${sample}.host_removed.fastq.gz ${ANALYSIS_OUTPUT}/${sample}/${sample}.clean.fastq.gz
       fi 
       # Find and link annotated VCF 
       if [ $(find ${ANALYSIS_OUTPUT}/${sample}/ -name "${sample}.pass.annotate.vcf" -type f) ] ; then 
            ln -sf ${sample}.pass.annotate.vcf ${ANALYSIS_OUTPUT}/${sample}/${sample}.pass.SnpEff.vcf
       fi 
    done
fi 


