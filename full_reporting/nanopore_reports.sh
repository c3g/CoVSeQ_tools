#!/bin/bash
set -eu -o pipefail

# LSPQ REPORTING PIPELINE
##################################################
# NANOPORE REPORTS PARSER
# Dumps all nanopore CoVSeQ reports into single files for easy parsing. 
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
echo "Preparing nanopore reports ${INPUT_NANOPORE_PATH} ..." 

####################
# Create output directory
mkdir -p ${OUTPUT_REPORT_PATH}

####################
# Internal metrics files 
## Define the metrics header
METRICS_FILE_HEADER="Run,analysis,sample,cons.perc.N,cons.len,cons.perc.GC,fq.size.pass,bam.perc.align,\
bam.mean.cov,bam.med.cov,bam.max.min.ratio,bam.perc.50x,bam.perc.100x,bam.perc.250x,bam.perc.500x,\
bam.perc.1000x,bam.perc.2000x"

####################
## Output internal metrics reports
echo "Preparing internal metrics..." 
## Overwrite previous report with header for new report
METRICS_FILE="${OUTPUT_REPORT_PATH}/internal_metrics.csv"
echo ${METRICS_FILE_HEADER} > ${METRICS_FILE}

## Search for all "legacy" internal metrics reports (pattern = report_metrics.csv)
grep -v "^sample," $(find ${INPUT_NANOPORE_PATH}/*/analysis/*/* -name "*metrics.csv" -type f) |\
   sed s:${INPUT_NANOPORE_PATH}/:: |\
   sed s:/analysis/:,:g |\
   sed s:/.*metrics.csv::g |\
   sed s/\:/,/g |\
   sed s:/:,: >> ${METRICS_FILE}

####################
## Prepare ncov_tools reports
echo "Preparing ncov_tools summary reports"
NCOV_TOOLS_SUMMARY_HEADER="Run\tanalysis\tsample\trun_name\tnum_consensus_snvs\tnum_consensus_n\
\tnum_consensus_iupac\tnum_variants_snvs\tnum_variants_indel\tnum_variants_indel_triplet\tmean_sequencing_depth\
\tmedian_sequencing_depth\tqpcr_ct\tcollection_date\tnum_weeks\tscaled_variants_snvs\tgenome_completeness\
\tqc_pass\tlineage\tlineage_notes\twatch_mutations\r"

NCOV_TOOLS_SUMMARY="${OUTPUT_REPORT_PATH}/ncov_tools_summary.tsv"
echo -e ${NCOV_TOOLS_SUMMARY_HEADER} > ${NCOV_TOOLS_SUMMARY}
grep -v "^sample" $(find ${INPUT_NANOPORE_PATH}/*/analysis/*/ncov_tools/ -name "*_summary_qc.tsv" -type f) |\
    sed s:${INPUT_NANOPORE_PATH}/:: |\
    sed s:/analysis/:'   ':g |\
    sed s:/ncov_tools/qc_reports/.*_summary_qc.tsv::g |\
    sed s/\:/'	'/ >> ${NCOV_TOOLS_SUMMARY} 

####################
## Prepare ncov_tools reports
echo "Preparing ncov_tools negative control reports." 
NCOV_TOOLS_NEG_CTRL_HEADER="Run\tanalysis\tfile\tqc\tgenome_covered_bases\tgenome_total_bases\
\tgenome_covered_fraction\tamplicons_detected"

NEG_CTRL_REPORT="${OUTPUT_REPORT_PATH}/neg_control_report.tsv"
echo -e ${NCOV_TOOLS_NEG_CTRL_HEADER} > ${NEG_CTRL_REPORT}
grep -v "^file" $(find ${INPUT_NANOPORE_PATH}/*/analysis/*/ncov_tools/ -name "*_negative_control_report.tsv" -type f) |\
    sed s:${INPUT_NANOPORE_PATH}/:: |\
    sed s:/analysis/:'   ':g |\
    sed s:/ncov_tools/qc_reports/.*_negative_control_report.tsv::g |\
    sed s/\:/'	'/ >> ${NEG_CTRL_REPORT} 
    
####################
## Prepare ncov_tools reports
echo "Preparing ncov_tools amplicon coverage reports." 
NCOV_TOOLS_AMPLICON_COV_HEADER="Run\tanalysis\tsample\tamp1\tamp2\tamp3\tamp4\tamp5\tamp6\tamp7\tamp8\tamp9\tamp10\tamp11\
\tamp12\tamp13\tamp14\tamp15\tamp16\tamp17\tamp18\tamp19\tamp20\tamp21\tamp22\tamp23\tamp24\tamp25\tamp26\tamp27\tamp28\
\tamp29\tamp30\tamp31\tamp32\tamp33\tamp34\tamp35\tamp36\tamp37\tamp38\tamp39\tamp40\tamp41\tamp42\tamp43\tamp44\tamp45\
\tamp46\tamp47\tamp48\tamp49\tamp50\tamp51\tamp52\tamp53\tamp54\tamp55\tamp56\tamp57\tamp58\tamp59\tamp60\tamp61\tamp62\
\tamp63\tamp64\tamp65\tamp66\tamp67\tamp68\tamp69\tamp70\tamp71\tamp72\tamp73\tamp74\tamp75\tamp76\tamp77\tamp78\tamp79\
\tamp80\tamp81\tamp82\tamp83\tamp84\tamp85\tamp86\tamp87\tamp88\tamp89\tamp90\tamp91\tamp92\tamp93\tamp94\tamp95\tamp96\
\tamp97\tamp98"

AMPLICON_COV_REPORT="${OUTPUT_REPORT_PATH}/amplicon_coverage_report.tsv"
echo -e ${NCOV_TOOLS_AMPLICON_COV_HEADER} > ${AMPLICON_COV_REPORT}
grep -v "^sample" $(find ${INPUT_NANOPORE_PATH}/*/analysis/*/ncov_tools/ -name "*_amplicon_coverage_table.tsv" -type f) |\
    sed s:${INPUT_NANOPORE_PATH}/:: |\
    sed s:/analysis/:'   ':g |\
    sed s:/ncov_tools/qc_analysis/.*_amplicon_coverage_table.tsv::g |\
    sed s/\:/'	'/ >> ${AMPLICON_COV_REPORT} 
    


