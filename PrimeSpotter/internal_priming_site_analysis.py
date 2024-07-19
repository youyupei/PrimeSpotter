import pysam
import pandas as pd

from internal_prim_filter import find_polyA_site, merge_gene_group
import parsers


class ProcessGTF:
    def __init__(self, gtf_file):
        self.gtf_data = self.parse_genecode_gtf(gtf_file)
    
    @classmethod
    def parse_genecode_gtf(gtf_file, feature_level="transcript"):
        """
        Parse the GTF file and return a DataFrame
        input:
        - gtf_file: str, path to the GTF file
        - feature_level: str, level to parse the GTF file (default: transcript)  
            transcript: consider only the gene and transcript level
            subtranscript: consider feature:"gene", "transcript", "CDS", "UTR", "exon" 
        """
        if feature_level not in ["transcript", "subtranscript"]:
            helper.error_msg(f"Invalid feature level: {feature_level}")
            sys.exit(1)
        
        gtf_data = pd.read_csv(gtf_file, sep="\t", comment="#", header=None).iloc[:, [0,2,3,4,6,8]]
        gtf_data.columns = ["Chromosome", "Feature", "Start", "End", "Strand", "Meta"]
        gtf_data['Gene_id'] = gtf_data.Meta.str.extract('gene_id "(.+?)";', expand=False)

        if feature_level == "transcript":
            gtf_data = gtf_data[gtf_data.Feature.isin(["gene", "transcript"])]

        elif feature_level == "subtranscript":
            gtf_data = gtf_data[gtf_data.Feature.isin(["gene", "transcript", "CDS", "UTR", "exon"])]
            gtf_data['Transcript_id'] = gtf_data.Meta.str.extract('transcri pt_id "(.+?)";', expand=False)

        del gtf_data['meta']

        return gtf_data
    
    @classmethod
    def get_exon_from_gene_group(grp):
        """
        Get the exon information from the gene group
        grp: GroupBy object from pandas, group by Gene_id
        """
        rst = []
        if not len(grp):
            return []
        for row in grp.itertuples():
            rst.append((row.Start, row.End))
            rst = list(set(rst))
        return rst
    

    @classmethod
    def process_gene(grp, genome_ref):
        """
        process the gene group (single process)
        """
        rst_stat = {}
        # read the fasta file
        fasta_file = pysam.FastaFile(genome_ref)
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

    args = parsers.parse_arg()

    # get information from gtf
    gtf_data = parse_genecode_gtf(args.gtf_file)

    # get exons
    exon_df = gtf_data[gtf_data.Feature=="exon"].groupby("Gene_id").apply(ProcessGTF.get_exon_from_gene_group)
    exon_df.name = 'exons'
    #
    # merge the 3' UTR information to the gtf_data
    # drop the transcript_id column
    gtf_data = gtf_data[gtf_data.Feature.isin(["gene", "transcript"])]
    gtf_data = gtf_data.drop("Transcript_id", axis=1)
    gtf_data = gtf_data.merge(exon_df, on="Gene_id", how="left")
    gtf_data['exons'] = gtf_data['exons'].fillna("").apply(list)

    poly_A_stats = gtf_data.groupby('Gene_id').apply(ProcessGTF.process_gene, genome_ref=args.genome_ref)
    poly_A_stats.to_csv("polyA_stats.csv")

