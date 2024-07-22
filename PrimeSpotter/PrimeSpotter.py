# python3 internal_prim_filter.py -b <bam_file> -g <genome_ref> -o <output_sam> -p <output_prefix> -s <output_summary> -t <tag> -f <gtf_file> -n <processes>
# take a bam file as input and output read coverage for each position in the genome

import sys
import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter

import helper
import parsers

class GtfGeneClassGenerator:
    """
    Take a pandas DataFrame and yield a GeneClass object per group
    """
    def __init__(self, gtf_file):
        self.gtf = gtf_file
        self.gtf_df = self.parse_genecode_gtf()
        self.iterator = self.gene_generator()
    def __iter__(self):
        self.iterator = self.gene_generator()  # Reset the generator
        return self
    
    def __next__(self):
        return next(self.iterator)

    def gene_generator(self):
        self.gene_grps = self.gtf_df.groupby('Gene_id')
        for index, group in self.gene_grps:
            gene_pd_serie = self.merge_gene_group(group)
            GeneClass._creation_locked = False
            gene = GeneClass(gene_pd_serie)
            yield gene
        GeneClass._creation_locked = True

    def merge_gene_group(self, gene_group):
        return self._merge_gene_group(gene_group)

    def parse_genecode_gtf(self):
        return self._parse_genecode_gtf(self.gtf)
    
    @staticmethod
    def _parse_genecode_gtf(gtf_file, feature_level="transcript"):
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

        del gtf_data['Meta']

        return gtf_data

    @staticmethod
    def _merge_gene_group(gene_group):
        """
        Merge a gene group into a single row with  TSS information as a list
        """
        # make sure the strand are identical for all rows, else raise an warning
        if gene_group.Strand.nunique() != 1:
            helper.warning_msg(f"Strand information is not consistent for gene {gene_group.Gene_id.iloc[0]}"
                            "This gene will be skipped.", printit = True)
            return None
        else:
            strand = gene_group.Strand.iloc[0]
            if strand not in ["+", "-"]:
                helper.warning_msg(f"Strand information is not valid for gene {gene_group.Gene_id.iloc[0]}"
                                "This gene will be skipped.", printit = True)
                return None
        # make sure only 1 row in the group has Feature == "gene"
        if sum(gene_group.Feature == "gene") != 1:
            helper.warning_msg(f"In GTF file, expected 1 gene with gene id {gene_group.Gene_id.iloc[0]}, but got {sum(gene_group.Feature == 'gene')}"
                            "This gene will be skipped.", printit = True)
            return None
        else:
            rst_series = gene_group[gene_group.Feature == "gene"].iloc[0]

            # add a new column to store the number of transcripts (concatenate the numbers separated by ",")
            if strand == "+":
                rst_series["TSS"] = [x for x in gene_group[gene_group.Feature == "transcript"].Start]
            else:
                rst_series["TSS"] = [x for x in gene_group[gene_group.Feature == "transcript"].End]
            return rst_series

class GeneClass:
    _creation_locked = True
    def __init__(self, series):
        """        
        pandas Series, gene information with columns:
            - Chromosome
            - Start
            - End
            - Strand
        """
        if not isinstance(series, pd.Series):
            raise ValueError("Expected a pandas Series object")
        self.series = series
        self.filtered_polyA_internal_binding_sites = None

    def __new__(cls, *args, **kwargs):
        """
        New instances can only be created by the GtfGeneClassGenerator
        """
        if getattr(cls, "_creation_locked", True):
            raise ValueError("Cannot create instances of this class directly")
        return super().__new__(cls)
    
    def __getattr__(self, attr):
        if attr in {'__getstate__', '__setstate__'}:
            return object.__getattr__(self, attr)
        return getattr(self.series, attr)
            
    
    def find_polyA_site(self, fasta_fn: str, inplace=True):
        """
        New attribute: filtered_polyA_internal_binding_sites will be added
        """
        all_polyA_site = self._find_all_polyA_site(self, fasta_fn)
        if inplace:
            self.filtered_polyA_internal_binding_sites = \
                [x for x in all_polyA_site if not any([x[0]<y<x[1] for y in self.TSS])]
        else:
            raise NotImplementedError("Not implemented yet")

    @staticmethod
    def _find_all_polyA_site(gene, fasta_fn: str, window_size=10, 
                    minAorT=8, merge_dist=3, flank_butter=25) -> None:
        """
        find the polyA site for the gene
        return a list of tuple, each tuple is the start and end of the polyA site
        ordered by the position in the gene (decending  for + strand, acending for - strand)
        """
        fasta_handle = pysam.FastaFile(fasta_fn)
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
        
        fasta_handle.close()
def process_gene(gene: GeneClass , 
                    bam_file_path: str, 
                    fasta_file_path: str,
                    sam_tag: str = "IP") -> tuple:

    # adding filtered_polyA_internal_binding_sites
    gene.find_polyA_site(fasta_file_path, inplace=True)
    
    # Initialise the return variables
    sam_buffer = ''
    counter = Counter()

    # read the bam
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")

    for read in bam_file.fetch(gene.Chromosome, gene.Start, gene.End):
        
        int_prim_flag = False
        if read.is_unmapped or read.is_supplementary:
            counter['skipped_reads'] += 1
            continue

        for start,end in gene.filtered_polyA_internal_binding_sites:
            if gene.Strand=='-':
                if start <= read.reference_start <= end:
                    int_prim_flag = True
                    
                    break
                elif read.reference_start < start:
                    break

            elif gene.Strand=='+':
                if start <= read.reference_end <= end:
                    int_prim_flag = True
                    counter['internal_priming'] += 1
                    break
                elif end < read.reference_end:
                    break
        

        if int_prim_flag:
            read.set_tag(sam_tag, 'Y', value_type='A')
            counter['internal_priming'] += 1
        else:
            read.set_tag(sam_tag, 'N', value_type='A')
            counter['normal_priming'] += 1

        sam_buffer+= read.to_string() + '\n'
    
    bam_file.close()

    return sam_buffer, counter

def write_buffer(sam_buffer, output_file=""):
    """
    Write the buffer to the output file
    """
    if output_file == "":
        # write to the stdout
        try:
            stdout.write(sam_buffer)
        except BrokenPipeError:
            sys.exit(0)
    else:
        with open(output_file, 'w') as f:
            f.write(sam_buffer)

def main():
    args = parsers.parse_arg()
    bam_template = pysam.AlignmentFile(args.bam_file, "rb")
    bam_header = bam_template.text
    write_buffer(bam_header, args.output_sam)

    # Expect counts in the counter:
        # - internal_priming
        # - normal_priming
        # - skipped_reads
    rst_counter = Counter()

    # Start
    gene_iter = GtfGeneClassGenerator(args.gtf_file)
    # Single process
    if args.processes == 1:
        for gene in tqdm(gene_iter):
            
            sam_buffer, counter = \
                process_gene(gene, args.bam_file, args.genome_ref)
            write_buffer(sam_buffer)
            rst_counter += counter

    # Multi process
    else:
        # release the lock for pickling in the multiprocessing
        # i=0
        # while i != None:
        #     next(gene_iter, None)
        # warning_msg("FINISH!", printit=True)
        for future in helper.multiprocessing_submit(
                            process_gene, 
                            gene_iter,
                            pbar=True,
                            pbar_unit = " gene_groups",
                            pbar_func = lambda x: 1,
                            n_process=args.processes, 
                            bam_file_path=args.bam_file, 
                            fasta_file_path=args.genome_ref):
        
            sam_buffer, counter = future.result()
            write_buffer(sam_buffer)
            rst_counter += counter
    
    # Finished
    print("Finished")


    with open(args.output_summary, 'w') as f:
        f.write(
            f"Total number of reads (mapped): {rst_counter['internal_priming'] + rst_counter['normal_priming']}\n"
            f"Internally primed reads: {rst_counter['internal_priming']}\n"
            f"Normally primed reads: {rst_counter['normal_priming']}\n"
            f"Proportion of internal priming reads: {rst_counter['internal_priming']/(rst_counter['internal_priming']+rst_counter['normal_priming'])}\n"
        )
        
if __name__ == "__main__":
    # redirect the print function to the stderror
    stdout = sys.stdout
    sys.stdout = sys.stderr
    main()