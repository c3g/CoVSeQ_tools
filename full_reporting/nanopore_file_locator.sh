#!/bin/bash
set -eu -o pipefail

# LSPQ REPORTING PIPELINE
##################################################
# NANOPORE FILE LOCATOR
# Locates output files for nanopore analyses
# Based on GenPipes `covseq` pipeline naming
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
echo "Searching for all nanopore output files in ${INPUT_NANOPORE_PATH} ..." 

####################
# Create output folders
mkdir -p ${OUTPUT_REPORT_PATH}/file_listings

####################
# Create run list
echo "Preparing run list..."

RUN_LIST=${OUTPUT_REPORT_PATH}/run_list.txt
echo -e "Run" > ${RUN_LIST}
find ${INPUT_NANOPORE_PATH} -mindepth 1 -maxdepth 1 -type d \
    -exec basename {} \; >> ${RUN_LIST}

####################
# Create dehosted fastq list
echo "Preparing dehosted fastq list..."

DEHOSTED_READ1=${OUTPUT_REPORT_PATH}/file_listings/dehosted_read1_list.txt
echo -e "filesize_dehosted_fastq_read1\tdehosted_fastq_read1" > ${DEHOSTED_READ1}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ \
   -name "*.clean.fastq*" -exec du -Dsh {} \; >> ${DEHOSTED_READ1}

####################
# Create alignment list
echo "Preparing alignment list..." 
ALIGNMENTS=${OUTPUT_REPORT_PATH}/file_listings/alignment_list.txt

echo -e "filesize_alignment_file\talignment_file" > ${ALIGNMENTS}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ -maxdepth 1 \
   -name "*.sorted.bam" -exec du -Dsh {} \; |\
   grep -v ".host_removed." |\
   grep -v ".hybrid." >> ${ALIGNMENTS}

####################
# Create nanoproe consensus list
echo "Preparing nanopolish consensus list..." 
NANOPOLISH_CONSENSUS=${OUTPUT_REPORT_PATH}/file_listings/nanopolish_consensus.txt

echo -e "filesize_nanopolish_consensus\tnanopolish_consensus" > ${NANOPOLISH_CONSENSUS}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ \
   -name "*.consensus.nanopore.*.fasta" -exec du -Dsh {} \; >> ${NANOPOLISH_CONSENSUS}

####################
# Create nanopolish VCF list
echo "Preparing nanopolish annotated VCF list..." 
NANOPOLISH_VCF=${OUTPUT_REPORT_PATH}/file_listings/nanopolish_vcf.txt

echo -e "filesize_nanopolish_annotated_vcf\tnanopolish_annotated_vcf" > ${NANOPOLISH_VCF}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ \
   -name "*.pass.SnpEff.vcf" -exec du -Dsh {} \; >> ${NANOPOLISH_VCF}

####################
echo "Finished searching for files in ${1}." 

