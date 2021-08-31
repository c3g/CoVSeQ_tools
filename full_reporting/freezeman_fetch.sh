#!/bin/bash
set -eu -o pipefail

# LSPQ REPORTING PIPELINE
##################################################
# FREEZEMAN FETCH SCRIPT
# Using the freezeman API, download the latest list of samples
# recorded on freezeman. Requires freezman username, a read-only 
# password file, and the CA BUNDLE permissions public key
##################################################

####################
if [ "$#" != 4 ]
then
    echo "Error: you provided an incorrect number of arguments." 
    echo "Usage: $0 <freezman_user> <password_file> <parth_to_cert> <output_dir>"
    exit 1
fi

####################
# Set file paths 
OUTPUT_DIR=${4}
mkdir -p ${OUTPUT_DIR}

####################
# Fetch latest freezeman report
echo "Fetching freezeman data..." 

python ${MUGQIC_INSTALL_HOME_DEV}/software/CoVSeQ_tools/CoVSeQ_tools-report_dev/full_reporting/freezeman_samples_export_by_cohort.py \
 -csv ${OUTPUT_DIR}/$(date +"%Y-%m-%d")_freezeman_fetch.csv \
 -cohort INSPQ_COVID \
 -url https://biobank.genome.mcgill.ca/api/ \
 -cert ${3} \
 -user ${1} -password $(cat ${2})

ln -sf $(date +"%Y-%m-%d")_freezeman_fetch.csv ${OUTPUT_DIR}/latest_freezeman_report.csv


