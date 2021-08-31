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
  echo "Usage: $0 <input_dir> <output_dir>"
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
# Create dehosted fastq list
echo "Preparing dehosted fastq list..."

DEHOSTED_READ1=${OUTPUT_REPORT_PATH}/file_listings/dehosted_read1_list.txt
echo "dehosted_fastq_read1" > ${DEHOSTED_READ1}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ \
   -type f -name "*.clean.fastq" >> ${DEHOSTED_READ1}

####################
# Create alignment list
echo "Preparing alignment list..." 
ALIGNMENTS=${OUTPUT_REPORT_PATH}/file_listings/alignment_list.txt

echo "alignment_file" > ${ALIGNMENTS}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ \
   -type f -name "*.trimmed.rg.sorted.bam" >> ${ALIGNMENTS}

####################
# Create nanoproe consensus list
echo "Preparing nanopolish consensus list..." 
NANOPOLISH_CONSENSUS=${OUTPUT_REPORT_PATH}/file_listings/nanopolish_consensus.txt

echo "nanopolish_consensus" > ${NANOPOLISH_CONSENSUS}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ \
   -type f -name "*.consensus.nanopore.*.fasta" >> ${NANOPOLISH_CONSENSUS}

####################
# Create nanopolish VCF list
echo "Preparing nanopolish annotated VCF list..." 
NANOPOLISH_VCF=${OUTPUT_REPORT_PATH}/file_listings/nanopolish_vcf.txt

echo "nanopolish_annotated_vcf" > ${NANOPOLISH_VCF}
find ${INPUT_NANOPORE_PATH}/*/analysis/*/*/ \
   -type f -name "*.pass.SnpEff.vcf" >> ${NANOPOLISH_VCF}

####################
echo "Finished searching for files in ${1}." 

