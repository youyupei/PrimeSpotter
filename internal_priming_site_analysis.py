import pysam
import numpy as np
import pandas as pd

# find plyA site
def find_polyA_site(gene, fasta_handle, window_size=10, 
                    minAorT=8, merge_dist=3, flank_butter=25):
    """
    find the polyA site for the gene
    return a list of tuple, each tuple is the start and end of the polyA site
    ordered by the position in the gene (decending  for + strand, acending for - strand)
    """
    seq = fasta_handle.fetch(gene.Chromosome, gene.Start, gene.End)
    if gene.Strand == "+":
        seq = np.array([int(x=='A') for x in seq])
    else:
        seq = np.array([int(x=='T') for x in seq])
    # get sliding window (size = 10) sum 
    window_sum = np.convolve(seq, np.ones(window_size), 'valid')
    # find location of consecutive window sum >= minAT
    polyA_site = np.where(window_sum >= minAorT)[0]
    if len(polyA_site):
        # merge consecutive polyA when the difference is <=3
        polyA_site = np.split(polyA_site, np.where(np.diff(polyA_site) > merge_dist)[0]+1)
        polyA_site = \
            [(gene.Start+x[0]-flank_butter, gene.Start+x[-1]+window_size+flank_butter) for x in polyA_site]
        if gene.Strand == "+":
            polyA_site = polyA_site[::-1]
        return polyA_site
    else:
        return []

def merge_gene_group(gene_group):
    # make sure the strand are identical for all rows, else raise an warning
    if gene_group.Strand.nunique() != 1:
        helper.warning_msg(f"Strand information is not consistent for gene {gene_group.Gene_id.iloc[0]}"
                           "This gene will be skipped.")
        return None
    else:
        strand = gene_group.Strand.iloc[0]
        if strand not in ["+", "-"]:
            helper.warning_msg(f"Strand information is not valid for gene {gene_group.Gene_id.iloc[0]}"
                               "This gene will be skipped.")
            return None
    # make sure only 1 row in the group has Feature == "gene"
    if sum(gene_group.Feature == "gene") != 1:
        helper.warning_msg(f"In GTF file, Expect 1 gene with gene id {gene_group.Gene_id.iloc[0]}, but got {sum(gene_group.Feature == 'gene')}"
                           "This gene will be skipped.")
        return None
    else:
        rst_df = gene_group[gene_group.Feature == "gene"]
        # deep copy the gene information
        rst_df = rst_df.copy(deep=True)
        # add a new column to store the number of transcripts (concatenate the numbers separated by ",")
        # 
        if strand == "+":
            rst_df["TTS"] = ','.join([str(x) for x in gene_group[gene_group.Feature == "transcript"].End])
        else:
            rst_df["TTS"] = ','.join([str(x) for x in gene_group[gene_group.Feature == "transcript"].Start])
        return rst_df

def get_exon_from_gene_group(grp):
    rst = []
    if not len(grp):
        return []
    for row in grp.itertuples():
        rst.append((row.Start, row.End))
        rst = list(set(rst))
    return rst


def parse_gtf(gtf_file):
    gtf_data = pd.read_csv(gtf_file, sep="\t", comment="#", header=None)
    #
    # remove unneeded information
    gtf_data = gtf_data[gtf_data[2].isin(["gene", "transcript", "CDS", "UTR", "exon"])]
    gtf_data = gtf_data.iloc[:, [0,2,3,4,6,8]]
    #
    # Extract the gene name from the last column
    gtf_data[9] = gtf_data[8].str.extract('transcript_id "(.+?)";', expand=False)
    gtf_data[8] = gtf_data[8].str.extract('gene_id "(.+?)";', expand=False)
    #
    # add column names
    gtf_data.columns = ["Chromosome", "Feature", "Start", "End", "Strand", "Gene_id", "transcript_id"]
    #
    # # get 3' UTR
    # utr3_df = gtf_data[gtf_data.Feature.isin(["CDS", "UTR"])].groupby(["Gene_id","transcript_id"])
    # utr3_df = utr3_df.apply(get_utr_from_transcript_group)
    # utr3_df = utr3_df.droplevel("transcript_id") 
    # utr3_df= utr3_df.groupby("Gene_id").apply(lambda y: sorted(list(set().union(*y))))
    # utr3_df.name = 'UTR3'
    #
    # get exons
    exon_df = gtf_data[gtf_data.Feature=="exon"].groupby("Gene_id").apply(get_exon_from_gene_group)
    exon_df.name = 'exons'
    #
    # merge the 3' UTR information to the gtf_data
    # drop the transcript_id column
    gtf_data = gtf_data[gtf_data.Feature.isin(["gene", "transcript"])]
    gtf_data = gtf_data.drop("transcript_id", axis=1)
    gtf_data = gtf_data.merge(exon_df, on="Gene_id", how="left")
    gtf_data['exons'] = gtf_data['exons'].fillna("").apply(list)
    #
    return gtf_data


def process_gene(grp, fasta_file_path):
    """
    process the gene group (single process)
    """
    rst_stat = {}
    # read the fasta file
    fasta_file = pysam.FastaFile(fasta_file_path)
    # get the gene information
    gene = merge_gene_group(grp).iloc[0]
    # get strand
    # reverse_strand = gene.Strand=='-'
    # find the polyA site
    polyA_site = find_polyA_site(gene, fasta_file)
    rst_stat['n_polyA_site'] = len(polyA_site) 
    # filter the polyA site that overlap with a list of location
    filter_loc = [int(x) for x in gene.TTS.split(",")]
    rst_stat['unique_TTS'] = len(filter_loc) 
    # # number of polyA site before filtering
    # num_polyA_site = len(polyA_site)
    # filter the polyA site that overlap with the Known TTS
    polyA_site = \
        [x for x in polyA_site if not any([x[0]<y<x[1] for y in filter_loc])]
    # get only the polyA site that overlap with the exon
    rst_stat['n_exonic_polyA_site'] = len([x for x in polyA_site if any(x[0]<y[1] and y[0]<x[1] for y in gene.exons)])
    return pd.Series(rst_stat)

if __name__ == "__main__":
    fasta_file_path = "/home/users/allstaff/you.yu/project/ref_files/GRCh38.primary_assembly.genome.fa"
    gtf_file = "/home/users/allstaff/you.yu/project/ref_files/gencode.v33.annotation.gtf"

    # get information from gtf
    gtf_data = parse_gtf(gtf_file)
    poly_A_stats = gtf_data.groupby('Gene_id').apply(process_gene, fasta_file_path=fasta_file_path)
    poly_A_stats.to_csv("polyA_stats.csv")

