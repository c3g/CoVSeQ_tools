#!/bin/bash
set -eu -o pipefail

# LSPQ REPORTING PIPELINE
##################################################
# ILLUMINA FILE LOCATOR
# Locates output files for illumina/MGI analyses
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
INPUT_ILLUMINA_PATH=${1}
OUTPUT_REPORT_PATH=${2}
echo "Searching for all illumina output files in ${INPUT_ILLUMINA_PATH} ..." 

####################
# Create output folders
mkdir -p ${OUTPUT_REPORT_PATH}/file_listings

####################
# Create dehosted fastq list
echo "Preparing dehosted fastq list..."

DEHOSTED_READ1=${OUTPUT_REPORT_PATH}/file_listings/dehosted_read1_list.txt
echo "dehosted_fastq_read1" > ${DEHOSTED_READ1}
find ${INPUT_ILLUMINA_PATH}/*/host_removal/* \
   -type f -name "*.host_removed.pair1.fastq.gz" >> ${DEHOSTED_READ1}
DEHOSTED_READ2=${OUTPUT_REPORT_PATH}/file_listings/dehosted_read2_list.txt
echo "dehosted_fastq_read1" > ${DEHOSTED_READ2}
find ${INPUT_ILLUMINA_PATH}/*/host_removal/* \
   -type f -name "*.host_removed.pair2.fastq.gz" >> ${DEHOSTED_READ2}

####################
# Create alignment list
echo "Preparing alignment list..." 
ALIGNMENTS=${OUTPUT_REPORT_PATH}/file_listings/alignment_list.txt

echo "alignment_file" > ${ALIGNMENTS}
find ${INPUT_ILLUMINA_PATH}/*/alignment/* \
   -type f -name "*.sorted.filtered.bam" >> ${ALIGNMENTS}

####################
# Create iVar consensus list
echo "Preparing iVar consensus list..." 
IVAR_CONSENSUS=${OUTPUT_REPORT_PATH}/file_listings/ivar_consensus.txt

echo "ivar_consensus" > ${IVAR_CONSENSUS}
find ${INPUT_ILLUMINA_PATH}/*/consensus/* \
   -type f -name "*.consensus.illumina.*.fasta" |\
   grep -v "freebayes_calling" >> ${IVAR_CONSENSUS}

####################
# Create FreeBayes consensus list
echo "Preparing FreeBayes consensus list..." 
FBAYES_CONSENSUS=${OUTPUT_REPORT_PATH}/file_listings/freebayes_consensus.txt

echo "freebayes_consensus" > ${FBAYES_CONSENSUS}
find ${INPUT_ILLUMINA_PATH}/*/consensus/* \
   -type f -name "*.freebayes_calling.consensus.illumina.*.fasta" \
   >> ${FBAYES_CONSENSUS}

####################
# Create iVar VCF list
echo "Preparing iVar annotated VCF list..." 
IVAR_VCF=${OUTPUT_REPORT_PATH}/file_listings/ivar_vcf.txt

echo "ivar_annotated_vcf" > ${IVAR_VCF}
find ${INPUT_ILLUMINA_PATH}/*/variant/* \
   -type f -name "*.sorted.filtered.primerTrim.annotate.vcf" >> ${IVAR_VCF}

####################
# Create FreeBayes VCF list
echo "Preparing FreeBayes annotated VCF list..." 
FBAYES_VCF=${OUTPUT_REPORT_PATH}/file_listings/freebayes_vcf.txt

echo "freebayes_annotated_vcf" > ${FBAYES_VCF}
find ${INPUT_ILLUMINA_PATH}/*/variant/* \
   -type f -name "*.freebayes_calling.fixed.norm.annotate.vcf" >> ${FBAYES_VCF}

####################
echo "Finished searching for files in ${1}." 

