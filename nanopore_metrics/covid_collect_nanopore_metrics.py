#!/usr/bin/env python
"""
Script that generates, colletcts and formats all relevant COVID metrics.

Arguments:
sample_name (str): the name of the sample used as input

Note: assumes that it is being run from the 'analysis' location of the sample, and that all files
follow the naming convention of our adapted ARTIC pipeline

Outputs a table with the follwing columns:
 - sample               = sample name
 - cons.perc.N          = percent N in consensus
 - cons.len             = length of consensus
 - cons.perc.GC         = GC content of consensus (%)
 - fq.size.pass         = Percent of fastq reads that pass size selection
 - bam.perc.align       = Percent of reads aligning to genome (from size pass)
 - bam.mean.cov         = Mean coverage
 - bam.med.cov          = Median coverage
 - bam.max.min.ratio    = Coverage max-min ratio
 - bam.perc.50x         = Percent covered at min 50x depth
 - bam.perc.100x        = Percent covered at min 100x depth
 - bam.perc.250x        = Percent covered at min 250x depth
 - bam.perc.500x        = Percent covered at min 500x depth
 - bam.perc.1000x       = Percent covered at min 1000x depth
 - bam.perc.2000x       = Percent covered at min 2000x depth
"""


import argparse
import pandas as pd
import pickle
import sys

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC

def parseoptions():
    """Command line options"""
    parser = argparse.ArgumentParser(description="Calculate the Jaccard coefficient within a sliding window from a bam alignment file")

    parser.add_argument('-c',
                        '--consensus',
                        help="Consensus file Fasta format.",
                        required=True)
    parser.add_argument('-s',
                        '--sample',
                        help="Sample name.",
                        required=True)
    parser.add_argument('-fqs',
                        '--fq_stats',
                        help="\"fastq.stats\" file.",
                        required=True)
    parser.add_argument('-pk',
                        '--pickle',
                        help="File format pickle.",
                        required=True)
    parser.add_argument('-o',
                        '--output',
                        help="Output filename.",
                        required=False)

    return parser.parse_args()

# Read arguments
sample = sys.argv[1]


def main():
    """main function"""

    # ARGS
    args = parseoptions()
    sample = args.sample
    ##########################################################################################
    # Construct file names
    consensus_file = args.consensus
    fastq_stats_file = args.fq_stats
    pickle_file = args.pickle
    output_file = args.output if args.output else None

    # Basic values
    CoV2_genome_size = float(29903)
    CoV2_chr_name = "MN908947.3"


    ##########################################################################################
    # Import Consensus sequence
    with open(consensus_file, 'r') as infile:
        consensus = SeqIO.read(infile, 'fasta')

    # Calculate consensus metrics
    cons_len = float(len(consensus.seq)) # Consensus length
    num_N = float(consensus.seq.count('N')) # Number of N in consensus
    perc_N = (num_N / cons_len) * 100 # Percent of N in consensus
    perc_GC = GC(consensus.seq)


    ##########################################################################################
    # Read in Fastq stats
    with open(fastq_stats_file, 'r') as infile:
        fastq_stats = pd.read_table(infile, sep=' ', header=None, index_col=0).transpose()


    # Calculate fastq metrics
    raw_reads = fastq_stats['raw_reads'][1]
    pass_reads = fastq_stats['pass_reads'][1]
    perc_pass = (pass_reads / raw_reads) * 100


    ##########################################################################################
    # Read in pickle file
    bam_pickle = pd.read_pickle(pickle_file)

    # Parse coverage info from pickle
    coverage_df = pd.DataFrame.from_dict(bam_pickle['pileup_stats']['coverage'][CoV2_chr_name], orient = 'index')

    # Calculate BAM metrics
    read_mapped = float(bam_pickle['read_stats']['mapped'])
    perc_mapped = (read_mapped / pass_reads) * 100
    mean_cov = coverage_df[0].mean()
    med_cov = coverage_df[0].median()
    max_cov = float(coverage_df[0].max())
    min_cov = float(coverage_df[0].min())
    max_min_ratio = (max_cov - min_cov) / mean_cov
    perc_50x = (len(coverage_df[coverage_df[0] > 50]) / CoV2_genome_size) * 100
    perc_100x = (len(coverage_df[coverage_df[0] > 100]) / CoV2_genome_size) * 100
    perc_250x = (len(coverage_df[coverage_df[0] > 250]) / CoV2_genome_size) * 100
    perc_500x = (len(coverage_df[coverage_df[0] > 500]) / CoV2_genome_size) * 100
    perc_1000x = (len(coverage_df[coverage_df[0] > 1000]) / CoV2_genome_size) * 100
    perc_2000x = (len(coverage_df[coverage_df[0] > 2000]) / CoV2_genome_size) * 100


    #######################################c
    # Create output dataframe

    output_dict = {
        "sample" : [sample],
        "cons.perc.N" : [perc_N],
        "cons.len" : [cons_len],
        "cons.perc.GC" : [perc_GC],
        "fq.size.pass" : [perc_pass],
        "bam.perc.align" : [perc_mapped],
        "bam.mean.cov" : [mean_cov],
        "bam.med.cov" : [med_cov],
        "bam.max.min.ratio" : [max_min_ratio],
        "bam.perc.50x" : [perc_50x],
        "bam.perc.100x" : [perc_100x],
        "bam.perc.250x" : [perc_250x],
        "bam.perc.500x" : [perc_500x],
        "bam.perc.1000x" : [perc_1000x],
        "bam.perc.2000x" : [perc_2000x]
    }

    output_df = pd.DataFrame.from_dict(output_dict, orient='columns')


    ##########################################################################################
    # Write output dataframe

    output_df.to_csv(output_file, sep = ',', index=False)

if __name__ == "__main__":
    main()
