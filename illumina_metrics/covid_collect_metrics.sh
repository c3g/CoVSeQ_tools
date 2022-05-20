#!/usr/bin/env bash

GENPIPES_OUTPUT_PATH=`pwd`

usage (){

echo
echo "usage: $0 <READSET_FILE>
  Gathers metrics for covseq pipeline from GenPipes. This script assumes you have samtools version 1.10 or above loaded in your environment"
echo
echo "   -t <THREADS>                  Number of threads to use (Default: 1)."
echo "   -r <READSET_FILE>             readset file used for GenPipes covseq analysis."
echo "   -o <GENPIPES_OUTPUT_PATH>     path of GenPipes covseq output location. (Default: $GENPIPES_OUTPUT_PATH)"

}

THREADS=1
while getopts "ht:r:o:" opt; do
  case $opt in
    t)
      THREADS=${OPTARG}
    ;;
    r)
      READSET_FILE=${OPTARG}
    ;;
    o)
      GENPIPES_OUTPUT_PATH=${OPTARG}
    ;;
    h)
      usage
      exit 0
    ;;
    \?)
      usage
      exit 1
    ;;
  esac
done


if [[ -z "$READSET_FILE" ]]; then
   usage
      exit 1
fi

METRICS_IVAR_OUT=${GENPIPES_OUTPUT_PATH}/metrics/metrics.csv
METRICS_FREEBAYES_OUT=${GENPIPES_OUTPUT_PATH}/metrics/metrics_freebayes.csv
HOST_CONTAMINATION_METRICS=${GENPIPES_OUTPUT_PATH}/metrics/host_contamination_metrics.tsv
HOST_REMOVED_METRICS=${GENPIPES_OUTPUT_PATH}/metrics/host_removed_metrics.tsv
KRAKEN_METRICS=${GENPIPES_OUTPUT_PATH}/metrics/kraken2_metrics.tsv

cat /dev/null > $METRICS_IVAR_OUT
cat /dev/null > $METRICS_FREEBAYES_OUT
cat /dev/null > $HOST_CONTAMINATION_METRICS
cat /dev/null > $HOST_REMOVED_METRICS
cat /dev/null > $KRAKEN_METRICS

echo "sample,cons.per.N,cons.len,cons.perc.GC,cons.perc.genome_frac,cons.N_per_100_kbp,fq.trim.pass,bam.perc.align,bam.filter.pass,bam.primertrim.pass,bam.mean.cov,bam.med.cov,bam.max-min/mean.cov,bam.perc.20x,bam.perc.50x,bam.perc.100x,bam.perc.250x,bam.perc.500x,bam.perc.1000x,bam.perc.2000x,bam.mean.insertsize,bam.med.insertsize,bam.sd.insertsize,bam.min.insertsize,bam.max.insertsize" > $METRICS_IVAR_OUT
echo "sample,cons.per.N,cons.len,cons.perc.GC,cons.perc.genome_frac,cons.N_per_100_kbp,fq.trim.pass,bam.perc.align,bam.filter.pass,bam.primertrim.pass,bam.mean.cov,bam.med.cov,bam.max-min/mean.cov,bam.perc.20x,bam.perc.50x,bam.perc.100x,bam.perc.250x,bam.perc.500x,bam.perc.1000x,bam.perc.2000x,bam.mean.insertsize,bam.med.insertsize,bam.sd.insertsize,bam.min.insertsize,bam.max.insertsize" > $METRICS_FREEBAYES_OUT

echo -e "Sample\tTotal_aligned\tHuman_only\tHuman_only_perc\tSARS_only\tSARS_only_perc\tUnmapped_only\tUnmapped_only_perc" > $HOST_CONTAMINATION_METRICS
echo -e "Sample\tTotal_aligned\tHuman_only\tHuman_only_perc\tSARS_only\tSARS_only_perc\tUnmapped_only\tUnmapped_only_perc" > $HOST_REMOVED_METRICS

echo -e "Sample\tHomo_sapiens_clade\tHomo_sapiens_clade_perc" > $KRAKEN_METRICS

for sample in `awk 'NR>1 {print $1}' $READSET_FILE`
do
    genome_size=29903
    consensus_fa=alignment/$sample/$sample.sorted.filtered.primerTrim.consensus.fa

    # Parsing quast results from ivar
    ivar_quast_tsv=`echo "metrics/dna/$sample/quast_metrics_ivar/report.tsv"`
    ivar_quast_html=`echo "metrics/dna/$sample/quast_metrics_ivar/report.html"`

    if [ -f "$ivar_quast_tsv" -a -f "$ivar_quast_html" ]; then
        ivar_cons_len=`grep -oP "Total length \(>= 0 bp\)\t\K.*?(?=$)" $ivar_quast_tsv`
        ivar_N_count=`grep -oP "# N's\",\"quality\":\"Less is better\",\"values\":\[\K.*?(?=])" $ivar_quast_html`
        if [ "$ivar_cons_len" != "0" ]; then
            ivar_cons_perc_N=`echo "scale=2; 100*$ivar_N_count/$ivar_cons_len" | bc -l`
        else
            ivar_cons_perc_N="NULL"
        fi
        ivar_cons_GC=`grep -oP "^GC \(%\)\t\K.*?(?=$)" $ivar_quast_tsv`
        ivar_cons_genome_frac=`grep -oP "^Genome fraction \(%\)\t\K.*?(?=$)" $ivar_quast_tsv`
        ivar_cons_N_perkbp=`grep -oP "^# N's per 100 kbp\t\K.*?(?=$)" $ivar_quast_tsv`
    else
        ivar_cons_len="NULL"
        ivar_N_count="NULL"
        ivar_cons_perc_N="NULL"
        ivar_cons_GC="NULL"
        ivar_cons_genome_frac="NULL"
        ivar_cons_N_perkbp="NULL"
    fi

    if [ -z $ivar_cons_GC ]; then
        ivar_cons_GC="NULL"
    fi
    if [ -z $ivar_cons_genome_frac ]; then
        ivar_cons_genome_frac="NULL"
    fi
    if [ -z $ivar_cons_N_perkbp ]; then
        ivar_cons_N_perkbp="NULL"
    fi

    # Parsing quast results from freebayes
    freebayes_quast_tsv=`echo "metrics/dna/$sample/quast_metrics_freebayes/report.tsv"`
    freebayes_quast_html=`echo "metrics/dna/$sample/quast_metrics_freebayes/report.html"`

    if [ -f "$freebayes_quast_tsv" -a -f "$freebayes_quast_html" ]; then
        freebayes_cons_len=`grep -oP "Total length \(>= 0 bp\)\t\K.*?(?=$)" $freebayes_quast_tsv`
        freebayes_N_count=`grep -oP "# N's\",\"quality\":\"Less is better\",\"values\":\[\K.*?(?=])" $freebayes_quast_html`
        if [ "$freebayes_cons_len" != "0" ]; then
            freebayes_cons_perc_N=`echo "scale=2; 100*$freebayes_N_count/$freebayes_cons_len" | bc -l`
        else
            freebayes_cons_perc_N="NULL"
        fi
        freebayes_cons_GC=`grep -oP "^GC \(%\)\t\K.*?(?=$)" $freebayes_quast_tsv`
        freebayes_cons_genome_frac=`grep -oP "^Genome fraction \(%\)\t\K.*?(?=$)" $freebayes_quast_tsv`
        freebayes_cons_N_perkbp=`grep -oP "^# N's per 100 kbp\t\K.*?(?=$)" $freebayes_quast_tsv`
    else
        freebayes_cons_len="NULL"
        freebayes_N_count="NULL"
        freebayes_cons_perc_N="NULL"
        freebayes_cons_GC="NULL"
        freebayes_cons_genome_frac="NULL"
        freebayes_cons_N_perkbp="NULL"
    fi

    if [ -z $freebayes_cons_GC ]; then
        freebayes_cons_GC="NULL"
    fi
    if [ -z $freebayes_cons_genome_frac ]; then
        freebayes_cons_genome_frac="NULL"
    fi
    if [ -z $freebayes_cons_N_perkbp ]; then
        freebayes_cons_N_perkbp="NULL"
    fi

    fq_surviving_trim=0
    tot_aligned_hybrid=0
    hum_only_hybrid=0
    sars_only_hybrid=0
    unmapped_only_hybrid=0
    tot_aligned_hrem=0
    hum_only_hrem=0
    sars_only_hrem=0
    unmapped_only_hrem=0
    homo_sapiens_clade=0
    homo_sapiens_clade_perc=0
    readset_count=0
    for readset_name in `grep "$sample" $READSET_FILE | awk '{print $2}'`
    do
        readset_count=$((readset_count+1))
        # Parsing cutadapt results
        cutadapt_file=`ls -t job_output/cutadapt/cutadapt.${readset_name}_*[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].o | head -n 1`
        fq_surviving_trim+=`grep -oP 'Pairs written \(passing filters\):.*\(\K.*?(?=%)' $cutadapt_file`

        # Computing host contamination and host cleaning metrics
        # samtools idxstats -@ $THREADS host_removal/${sample}/${readset_name}*.hybrid.sorted.bam | awk -v sample=$sample '{array[$1]=$3; pwet[$1]=$4; next} END {tot=0; unmapped=0; for (chr in array) {tot+=array[chr]; unmapped+=pwet[chr]}; tot+=unmapped; cov=array["MN908947.3"]; hum=tot-cov-unmapped; if (tot <= 0) {printf sample"\t"tot"\t"hum"\tNULL\t"cov"\tNULL\t"unmapped"\tNULL\n"} else {printf sample"\t"tot"\t"hum"\t%.2f\t"cov"\t%.2f\t"unmapped"\t%.2f\n", 100*hum/tot, 100*cov/tot, 100*unmapped/tot}}' >> $HOST_CONTAMINATION_METRICS
        hybrid=`samtools idxstats -@ $THREADS host_removal/${sample}/${readset_name}*.hybrid.sorted.bam | awk '{array[$1]=$3; pwet[$1]=$4; next} END {tot=0; unmapped=0; for (chr in array) {tot+=array[chr]; unmapped+=pwet[chr]}; tot+=unmapped; cov=array["MN908947.3"]; hum=tot-cov-unmapped; printf tot","hum","cov","unmapped}'`
        tot_aligned_hybrid=$((tot_aligned_hybrid+$(echo $hybrid | cut -d$',' -f1)))
        hum_only_hybrid=$((hum_only_hybrid+$(echo $hybrid | cut -d$',' -f2)))
        sars_only_hybrid=$((sars_only_hybrid+$(echo $hybrid | cut -d$',' -f3)))
        unmapped_only_hybrid=$((unmapped_only_hybrid+$(echo $hybrid | cut -d$',' -f4)))

        hrem=`samtools idxstats -@ $THREADS host_removal/${sample}/${readset_name}*.host_removed.sorted.bam | awk '{array[$1]=$3; pwet[$1]=$4; next} END {tot=0; unmapped=0; for (chr in array) {tot+=array[chr]; unmapped+=pwet[chr]}; tot+=unmapped; cov=array["MN908947.3"]; hum=tot-cov-unmapped; printf tot","hum","cov","unmapped}'`
        tot_aligned_hrem=$((tot_aligned_hrem+$(echo $hybrid | cut -d$',' -f1)))
        hum_only_hrem=$((hum_only_hrem+$(echo $hybrid | cut -d$',' -f2)))
        sars_only_hrem=$((sars_only_hrem+$(echo $hybrid | cut -d$',' -f3)))
        unmapped_only_hrem=$((unmapped_only_hrem+$(echo $hybrid | cut -d$',' -f4)))
        # samtools idxstats -@ $THREADS host_removal/${sample}/${readset_name}*.host_removed.sorted.bam | awk -v sample=$sample '{array[$1]=$3; pwet[$1]=$4; next} END {tot=0; unmapped=0; for (chr in array) {tot+=array[chr]; unmapped+=pwet[chr]}; tot+=unmapped; cov=array["MN908947.3"]; hum=tot-cov-unmapped; if (tot <= 0) {printf sample"\t"tot"\t"hum"\tNULL\t"cov"\tNULL\t"unmapped"\tNULL\n"} else {printf sample"\t"tot"\t"hum"\t%.2f\t"cov"\t%.2f\t"unmapped"\t%.2f\n", 100*hum/tot, 100*cov/tot, 100*unmapped/tot}}' >> $HOST_REMOVED_METRICS

        kraken_file=`ls metrics/dna/${sample}/kraken_metrics/${readset_name}*.kraken2_report`
        if [ -s "$kraken_file" ]; then
            homo_sapiens_kraken=`grep "Homo sapiens" $kraken_file`
            homo_sapiens_clade=$((homo_sapiens_clade+$(echo $homo_sapiens_kraken | cut -d$',' -f2)))
            homo_sapiens_clade_perc=$((homo_sapiens_clade_perc+$(echo $homo_sapiens_kraken | cut -d$',' -f1)))
            # grep "Homo sapiens" $kraken_file | awk -v sample=$sample '{print sample"\t"$2"\t"$1}' >> $KRAKEN_METRICS

        else
            echo -e "$sample\tNULL\tNULL" >> $KRAKEN_METRICS
        fi
    done

    if ((tot_aligned_hybrid <= 0)); then
        printf $sample"\t"$tot_aligned_hybrid"\t"$hum_only_hybrid"\tNULL\t"$sars_only_hybrid"\tNULL\t"$unmapped_only_hybrid"\tNULL\n" >> $HOST_CONTAMINATION_METRICS
    else
        printf $sample"\t"$tot_aligned_hybrid"\t"$hum_only_hybrid"\t%.2f\t"$sars_only_hybrid"\t%.2f\t"$unmapped_only_hybrid"\t%.2f\n", $(echo "100*$hum_only_hybrid/$tot_aligned_hybrid" | bc -l), $(echo "100*$sars_only_hybrid/$tot_aligned_hybrid" | bc -l), $(echo "100*$unmapped_only_hybrid/$tot_aligned_hybrid" | bc -l) >> $HOST_CONTAMINATION_METRICS
    fi

    if ((tot_aligned_hrem <= 0)); then
        printf $sample"\t"$tot_aligned_hrem"\t"$hum_only_hrem"\tNULL\t"$sars_only_hrem"\tNULL\t"$unmapped_only_hrem"\tNULL\n" >> $HOST_REMOVED_METRICS
    else
        printf $sample"\t"$tot_aligned_hrem"\t"$hum_only_hrem"\t%.2f\t"$sars_only_hrem"\t%.2f\t"$unmapped_only_hrem"\t%.2f\n", $(echo "100*$hum_only_hrem/$tot_aligned_hrem" | bc -l), $(echo "100*$sars_only_hrem/$tot_aligned_hrem" | bc -l), $(echo "100*$unmapped_only_hrem/$tot_aligned_hrem" | bc -l) >> $HOST_REMOVED_METRICS
    fi

    printf $sample"\t"$homo_sapiens_clade"\t%.2f", $(echo "$homo_sapiens_clade_perc/$readset_count" | bc -l) >> $KRAKEN_METRICS

    # Parsing sambamba flagstat and picard metrics
    raw_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.flagstat"`
    filtered_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.filtered.flagstat"`
    primertrim_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.filtered.primerTrim.flagstat"`
    total_raw=`head -n 1 $raw_flagstat_file | grep -oP '^\K.*?(?= )'`
    alignment_metrics_file=`echo "alignment/$sample/$sample.sorted.filtered.all.metrics.alignment_summary_metrics"`
    total_filter=`head -n 1 $filtered_flagstat_file | grep -oP '^\K.*?(?= )'`
    if [ "$total_raw" != "0" -a "$total_filter" != "0" -a -f "$alignment_metrics_file" ]; then
        bam_surviving_filter=`echo "scale=2; 100*$total_filter/$total_raw" | bc -l`
        total_primertrim=`head -n 1 $primertrim_flagstat_file | grep -oP '^\K.*?(?= )'`
        bam_surviving_primertrim=`echo "scale=2; 100*$total_primertrim/$total_filter" | bc -l`
    else
        bam_surviving_filter="NULL"
        bam_surviving_primertrim="NULL"
    fi
    bam_aln=`grep "mapped (" $raw_flagstat_file | grep -oP '^.*\(\K.*?(?=%)'`

    # Parsing BedGraph file
    bam_bedgraph_file=`echo "alignment/$sample/$sample.sorted.filtered.BedGraph"`
    bam_meancov=`awk '{for (i = 1; i <= $3-$2; ++i) sum+=$4} END {print sum}' $bam_bedgraph_file`
    bam_meancov=`echo "scale=2; $bam_meancov/$genome_size" | bc -l`
    bam_mediancov=`awk '{for (i = 1; i <= $3-$2; ++i)  print $4}' $bam_bedgraph_file | sort -g | awk '{count[NR] = $1} END {if (NR % 2) {print count[(NR + 1) / 2];} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;}}'`
    bam_mincov=`awk '{print $4}' $bam_bedgraph_file | sort -g | head -n 1 | awk '{printf "%i", $0}'`
    bam_maxcov=`awk '{print $4}' $bam_bedgraph_file | sort -g | tail -n 1 | awk '{printf "%i", $0}'`
    if [ "$bam_meancov" == "0" ]; then
        bam_maxmincovmean="NULL"
    else
        bam_maxmincovmean=`echo "scale=2; ($bam_maxcov-$bam_mincov)/$bam_meancov" | bc -l`
    fi
    bam_cov20X=`awk '{if ($4 > 20) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov20X=`echo "scale=2; 100*$bam_cov20X/$genome_size" | bc -l`
    bam_cov50X=`awk '{if ($4 > 50) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
    bam_cov100X=`awk '{if ($4 > 100) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov100X=`echo "scale=2; 100*$bam_cov100X/$genome_size" | bc -l`
    bam_cov250X=`awk '{if ($4 > 250) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov250X=`echo "scale=2; 100*$bam_cov250X/$genome_size" | bc -l`
    bam_cov500X=`awk '{if ($4 > 500) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov500X=`echo "scale=2; 100*$bam_cov500X/$genome_size" | bc -l`
    bam_cov1000X=`awk '{if ($4 > 1000) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov1000X=`echo "scale=2; 100*$bam_cov1000X/$genome_size" | bc -l`
    bam_cov2000X=`awk '{if ($4 > 2000) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov2000X=`echo "scale=2; 100*$bam_cov2000X/$genome_size" | bc -l`

    # Parsing Picard insert size metrics
    bam_insertsize_file=`echo "metrics/dna/$sample/picard_metrics/$sample.sorted.all.metrics.insert_size_metrics"`
    if [ -f "$bam_insertsize_file" ]; then
        bam_meaninsertsize=`awk 'NR==8 {print $6}' $bam_insertsize_file`
        bam_meaninsertsize=`echo "scale=2; ${bam_meaninsertsize//,}/1" | bc -l`
        bam_medianinsertsize=`awk 'NR==8 {print $1}' $bam_insertsize_file`
        bam_mininsertsize=`awk 'NR==8 {print $4}' $bam_insertsize_file`
        bam_maxinsertsize=`awk 'NR==8 {print $5}' $bam_insertsize_file`
        bam_sdinsertsize=`awk 'NR==8 {print $7}' $bam_insertsize_file`
        if [ "$bam_sdinsertsize" != "?" ]; then
            bam_sdinsertsize=`echo "scale=2; ${bam_sdinsertsize//,}/1" | bc -l`
        else
            bam_sdinsertsize="NULL"
        fi
    else 
        bam_meaninsertsize="NULL"
        bam_medianinsertsize="NULL"
        bam_mininsertsize="NULL"
        bam_maxinsertsize="NULL"
        bam_sdinsertsize="NULL"
    fi

    echo "$sample,$ivar_cons_perc_N,$ivar_cons_len,$ivar_cons_GC,$ivar_cons_genome_frac,$ivar_cons_N_perkbp,$fq_surviving_trim,$bam_aln,$bam_surviving_filter,$bam_surviving_primertrim,$bam_meancov,$bam_mediancov,$bam_maxmincovmean,$bam_cov20X,$bam_cov50X,$bam_cov100X,$bam_cov250X,$bam_cov500X,$bam_cov1000X,$bam_cov2000X,$bam_meaninsertsize,$bam_medianinsertsize,$bam_sdinsertsize,$bam_mininsertsize,$bam_maxinsertsize" >> $METRICS_IVAR_OUT
    echo "$sample,$freebayes_cons_perc_N,$freebayes_cons_len,$freebayes_cons_GC,$freebayes_cons_genome_frac,$freebayes_cons_N_perkbp,$fq_surviving_trim,$bam_aln,$bam_surviving_filter,$bam_surviving_primertrim,$bam_meancov,$bam_mediancov,$bam_maxmincovmean,$bam_cov20X,$bam_cov50X,$bam_cov100X,$bam_cov250X,$bam_cov500X,$bam_cov1000X,$bam_cov2000X,$bam_meaninsertsize,$bam_medianinsertsize,$bam_sdinsertsize,$bam_mininsertsize,$bam_maxinsertsize" >> $METRICS_FREEBAYES_OUT
    # break
done
