#!/bin/bash
set -eu -o pipefail

# LSPQ REPORTING PIPELINE
##################################################
# ILLUMINA REPORTS PARSER
# Dumps all illumina CoVSeQ reports into single files for easy parsing. 
# Assumes data was analyzed with GenPipes `covseq` pipeline
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
echo "Preparing illumina reports ${INPUT_ILLUMINA_PATH} ..." 

####################
# Create output directory
mkdir -p ${OUTPUT_REPORT_PATH}

####################
# Internal metrics files 
## Define the metrics header
METRICS_FILE_HEADER="Run,Sample,Nb reads,Percent human reads,Nb clean reads,Mean coverage,Percent N,\
Length low cov region (<20X),Percent consensus > 100X,Length consensus,Nb variants > 10 perc allele freq,\
Nb variants > 75 perc allele freq,PASS/FLAG/REJ"

####################
## Output iVar internal metrics reports
echo "Preparing iVar internal metrics..." 
## Overwrite previous report with header for new report
IVAR_METRICS_FILE="${OUTPUT_REPORT_PATH}/ivar_internal_metrics.csv"
echo ${METRICS_FILE_HEADER} > ${IVAR_METRICS_FILE}

## Search for all "legacy" iVar internal metrics reports (pattern = report_metrics.csv)
grep -v "^Sample," $(find ${INPUT_ILLUMINA_PATH}/*/report -name "report_metrics.csv" -type f) |\
   sed s:${INPUT_ILLUMINA_PATH}/:: |\
   sed s:/report/report_metrics.csv::g |\
   sed s/\:/,/g >> ${IVAR_METRICS_FILE}
## Search for all "new" iVar internal metrics reprots (pattern = report_metrics_ivar.csv)
grep -v "^Sample," $(find ${INPUT_ILLUMINA_PATH}/*/report -name "report_metrics_ivar.csv" -type f) |\
   sed s:${INPUT_ILLUMINA_PATH}/:: |\
   sed s:/report/report_metrics_ivar.csv::g |\
   sed s/\:/,/g >> ${IVAR_METRICS_FILE}

####################
## Output Freebayes internal metrics reports
echo "Preparing Freebayes internal metrics" 
## Overwrite previous report with header for new report
FBAYES_METRICS_FILE="${OUTPUT_REPORT_PATH}/fbayes_internal_metrics.csv"
echo ${METRICS_FILE_HEADER} > ${FBAYES_METRICS_FILE}

## Search for all "new" internal metrics reprots (pattern = report_metrics_freebayes.csv)
grep -v "^Sample," $(find ${INPUT_ILLUMINA_PATH}/*/report -name "report_metrics_freebayes.csv" -type f) |\
   sed s:${INPUT_ILLUMINA_PATH}/:: |\
   sed s:/report/report_metrics_freebayes.csv::g |\
   sed s/\:/,/g >> ${FBAYES_METRICS_FILE}

####################
## Prepare ncov_tools reports
echo "Preparing ncov_tools summary reports"
NCOV_TOOLS_SUMMARY_HEADER="Run\tsample\trun_name\tnum_consensus_snvs\tnum_consensus_n\tnum_consensus_iupac\
\tnum_variants_snvs\tnum_variants_indel\tnum_variants_indel_triplet\tmean_sequencing_depth\
\tmedian_sequencing_depth\tqpcr_ct\tcollection_date\tnum_weeks\tscaled_variants_snvs\tgenome_completeness\tqc_pass\tlineage\tlineage_notes\twatch_mutations\r"

IVAR_NCOV_TOOLS_SUMMARY="${OUTPUT_REPORT_PATH}/ivar_ncov_tools_summary.tsv"
echo -e ${NCOV_TOOLS_SUMMARY_HEADER} > ${IVAR_NCOV_TOOLS_SUMMARY}
grep -v "^sample" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools/ -name "*_summary_qc.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools/qc_reports/.*_summary_qc.tsv::g |\
    sed s/\:/'	'/ >> ${IVAR_NCOV_TOOLS_SUMMARY} 
grep -v "^sample" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools_ivar/ -name "*_summary_qc.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools_ivar/qc_reports/.*_summary_qc.tsv::g |\
    sed s/\:/'	'/ >> ${IVAR_NCOV_TOOLS_SUMMARY} 

FBAYES_NCOV_TOOLS_SUMMARY="${OUTPUT_REPORT_PATH}/fbayes_ncov_tools_summary.tsv"
echo -e ${NCOV_TOOLS_SUMMARY_HEADER} > ${FBAYES_NCOV_TOOLS_SUMMARY}
grep -v "^sample" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools_freebayes/ -name "*_summary_qc.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools_freebayes/qc_reports/.*_summary_qc.tsv::g |\
    sed s/\:/'	'/ >> ${FBAYES_NCOV_TOOLS_SUMMARY} 

####################
## Prepare ncov_tools reports
echo "Preparing ncov_tools negative control reports." 
NCOV_TOOLS_NEG_CTRL_HEADER="Run\tfile\tqc\tgenome_covered_bases\tgenome_total_bases\tgenome_covered_fraction\tamplicons_detected"

IVAR_NEG_CTRL_REPORT="${OUTPUT_REPORT_PATH}/ivar_neg_control_report.tsv"
echo -e ${NCOV_TOOLS_NEG_CTRL_HEADER} > ${IVAR_NEG_CTRL_REPORT}
grep -v "^file" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools/ -name "*_negative_control_report.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools/qc_reports/.*_negative_control_report.tsv::g |\
    sed s/\:/'	'/ >> ${IVAR_NEG_CTRL_REPORT} 
grep -v "^file" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools_ivar/ -name "*_negative_control_report.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools_ivar/qc_reports/.*_negative_control_report.tsv::g |\
    sed s/\:/'	'/ >> ${IVAR_NEG_CTRL_REPORT} 
    
FBAYES_NEG_CTRL_REPORT="${OUTPUT_REPORT_PATH}/freebayes_neg_control_report.tsv"
echo -e ${NCOV_TOOLS_NEG_CTRL_HEADER} > ${FBAYES_NEG_CTRL_REPORT}
grep -v "^file" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools_freebayes/ -name "*_negative_control_report.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools_freebayes/qc_reports/.*_negative_control_report.tsv::g |\
    sed s/\:/'	'/ >> ${FBAYES_NEG_CTRL_REPORT} 

####################
## Prepare ncov_tools reports
echo "Preparing ncov_tools amplicon coverage reports." 
NCOV_TOOLS_AMPLICON_COV_HEADER="Run\tsample\tamp1\tamp2\tamp3\tamp4\tamp5\tamp6\tamp7\tamp8\tamp9\tamp10\tamp11\tamp12\t\
amp13\tamp14\tamp15\tamp16\tamp17\tamp18\tamp19\tamp20\tamp21\tamp22\tamp23\tamp24\tamp25\tamp26\tamp27\tamp28\tamp29\t\
amp30\tamp31\tamp32\tamp33\tamp34\tamp35\tamp36\tamp37\tamp38\tamp39\tamp40\tamp41\tamp42\tamp43\tamp44\tamp45\tamp46\t\
amp47\tamp48\tamp49\tamp50\tamp51\tamp52\tamp53\tamp54\tamp55\tamp56\tamp57\tamp58\tamp59\tamp60\tamp61\tamp62\tamp63\t\
amp64\tamp65\tamp66\tamp67\tamp68\tamp69\tamp70\tamp71\tamp72\tamp73\tamp74\tamp75\tamp76\tamp77\tamp78\tamp79\tamp80\t\
amp81\tamp82\tamp83\tamp84\tamp85\tamp86\tamp87\tamp88\tamp89\tamp90\tamp91\tamp92\tamp93\tamp94\tamp95\tamp96\tamp97\tamp98"

IVAR_AMPLICON_COV_REPORT="${OUTPUT_REPORT_PATH}/ivar_amplicon_coverage_report.tsv"
echo -e ${NCOV_TOOLS_AMPLICON_COV_HEADER} > ${IVAR_AMPLICON_COV_REPORT}
grep -v "^sample" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools/ -name "*_amplicon_coverage_table.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools/qc_analysis/.*_amplicon_coverage_table.tsv::g |\
    sed s/\:/'	'/ >> ${IVAR_AMPLICON_COV_REPORT} 
grep -v "^sample" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools_ivar/ -name "*_amplicon_coverage_table.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools_ivar/qc_analysis/.*_amplicon_coverage_table.tsv::g |\
    sed s/\:/'	'/ >> ${IVAR_AMPLICON_COV_REPORT} 
    
FBAYES_AMPLICON_COV_REPORT="${OUTPUT_REPORT_PATH}/freebayes_amplicon_coverage_report.tsv"
echo -e ${NCOV_TOOLS_AMPLICON_COV_HEADER} > ${FBAYES_AMPLICON_COV_REPORT}
grep -v "^sample" $(find ${INPUT_ILLUMINA_PATH}/*/report/ncov_tools_freebayes/ -name "*_amplicon_coverage_table.tsv" -type f) |\
    sed s:${INPUT_ILLUMINA_PATH}/:: |\
    sed s:/report/ncov_tools_freebayes/qc_analysis/.*_amplicon_coverage_table.tsv::g |\
    sed s/\:/'	'/ >> ${FBAYES_AMPLICON_COV_REPORT} 


