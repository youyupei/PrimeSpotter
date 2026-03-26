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
    def __init__(self, gtf_file, gene_subset=None):
        """
        Parameters
        ----------
        gtf_file : str
            Path to the GTF/GFF annotation file.
        gene_subset : set or None
            Optional set of Gene_ids to restrict processing to.
            If None, all genes in the GTF are processed.
        """
        self.gtf = gtf_file
        self.gene_subset = gene_subset
        self.gtf_df = self.parse_genecode_gtf() # the coordinates in the gtf file will be converted to 0-based
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
            # merge_gene_group returns None for malformed GTF rows — skip them
            if gene_pd_serie is None:
                continue
            GeneClass._creation_locked = False
            gene = GeneClass(gene_pd_serie)
            yield gene
        GeneClass._creation_locked = True

    def merge_gene_group(self, gene_group):
        return self._merge_gene_group(gene_group)

    def parse_genecode_gtf(self):
        gtf_data = self._parse_genecode_gtf(self.gtf)
        if self.gene_subset is not None:
            n_before = gtf_data['Gene_id'].nunique()
            # Match against gene name first, fall back to gene ID
            # This lets the gene list contain either KCNQ1OT1 or ENSG00000269821
            matched_by_name = gtf_data['Gene_name'].isin(self.gene_subset)
            matched_by_id   = gtf_data['Gene_id'].isin(self.gene_subset)
            gtf_data = gtf_data[matched_by_name | matched_by_id]
            n_after = gtf_data['Gene_id'].nunique()
            helper.green_msg(
                f"Gene subset applied: {n_after}/{n_before} genes retained.",
                printit=True
            )
            # Report anything from the list that matched neither name nor id
            matched_names = set(gtf_data['Gene_name'].dropna().unique())
            matched_ids   = set(gtf_data['Gene_id'].dropna().unique())
            missing = self.gene_subset - matched_names - matched_ids
            if missing:
                helper.warning_msg(
                    f"{len(missing)} entry/entries from --gene-list were not found in the GTF "
                    f"(checked both gene_name and gene_id): "
                    f"{', '.join(sorted(missing)[:10])}"
                    f"{'...' if len(missing) > 10 else ''}",
                    printit=True
                )
        return gtf_data
    
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
        # update 1-based to 0-based coordinate for start position
        gtf_data['Start'] = gtf_data['Start'] - 1
        gtf_data['Gene_id'] = gtf_data.Meta.str.extract('gene_id "(.+?)";', expand=False)
        gtf_data['Gene_name'] = gtf_data.Meta.str.extract('gene_name "(.+?)";', expand=False)

        if feature_level == "transcript":
            gtf_data = gtf_data[gtf_data.Feature.isin(["gene", "transcript"])]

        elif feature_level == "subtranscript":
            gtf_data = gtf_data[gtf_data.Feature.isin(["gene", "transcript", "CDS", "UTR", "exon"])]
            gtf_data['Transcript_id'] = gtf_data.Meta.str.extract('transcript_id "(.+?)";', expand=False)
            # update 1-based to 0-based coordinate for start position
            gtf_data['Start'] = gtf_data['Start'] - 1

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
            if strand == "+":
                rst_series["TSS"] = [x for x in gene_group[gene_group.Feature == "transcript"].Start]
                rst_series["TTS"] = [x for x in gene_group[gene_group.Feature == "transcript"].End]
            elif strand == "-":
                rst_series["TSS"] = [x for x in gene_group[gene_group.Feature == "transcript"].End]
                rst_series["TTS"] = [x for x in gene_group[gene_group.Feature == "transcript"].Start]
            else:
                raise ValueError(f"Invalid strand {strand} for gene {rst_series.Gene_id}")

            return rst_series

class GeneClass:
    _creation_locked = True
    def __init__(self, series):
        """        
        pandas Series, gene information with columns:
            - Chromosome
            - Start (0-based)
            - End (0-based)
            - Strand
            - Gene_id
            - TSS (list of transcript start sites) (0-based)
            - TTS (list of transcript termination sites) (0-based)
        """
        if not isinstance(series, pd.Series):
            raise ValueError("Expected a pandas Series object")
        self.series = series
        self.polyA_internal_binding_sites = None
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
            
    
    def find_polyA_site(self, fasta_fn: str, inplace=True, flank_buffer: int = 10):
        """
        Find the polyA site for the gene and store it in the filtered_polyA_internal_binding_sites attribute.
        """
        all_polyA_site = self._find_all_polyA_site(self, fasta_fn, flank_buffer=flank_buffer)
        if inplace:
            self.polyA_internal_binding_sites = all_polyA_site
            self.filtered_polyA_internal_binding_sites = \
                [x for x in all_polyA_site if not any([x[0]<y<x[1] for y in self.TTS])]
        else:
            raise NotImplementedError("Not implemented yet")

    @staticmethod
    def _find_all_polyA_site(gene, fasta_fn: str, window_size=10, 
                    minAorT=8, merge_dist=2, flank_buffer=10) -> None:
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
            # merge polyA windows with a when the gap between windows is <= merge_dist
            # note that polyA_site is the index of the start of the window so the correct calculation of the gap: site distance - window_size 
            polyA_site = np.split(polyA_site, np.where(np.diff(polyA_site) - window_size> merge_dist)[0]+1)
            polyA_site = \
                [(gene.Start+x[0]-flank_buffer, gene.Start+x[-1]+window_size+flank_buffer) for x in polyA_site]
            
            # forgot why I reversed the polyA_site for + strand, but it does not matter
            if gene.Strand == "+":
                polyA_site = polyA_site[::-1]

            return polyA_site
        else:
            return []
        
        fasta_handle.close()
def process_gene(gene: GeneClass , 
                    bam_file_path: str, 
                    fasta_file_path: str,
                    flank_buffer: int = 10) -> tuple:

    """
    Process the gene and return the offset of the reads in each category
    Input:
        - gene: GeneClass object
        - bam_file_path: str, path to the bam file
        - fasta_file_path: str, path to the fasta file
        - flank_buffer: int, bases to expand each poly-A site on both sides (default 10)
    Output:
        - int_prim_offsets: list of offsets for internal priming
        - normal_prim_offsets: list of offsets for normal priming
        - skipped_reads_offsets: list of offsets for skipped reads
    Note:
        The reason for output the offsets as list is that same read can reappear when processing
        different genes, so we need to store the offsets for each gene separately and merge them later
    """
    # adding filtered_polyA_internal_binding_sites and polyA_internal_binding_sites
    gene.find_polyA_site(fasta_file_path, inplace=True, flank_buffer=flank_buffer)
    n_site_filtered = len(gene.filtered_polyA_internal_binding_sites)
    n_site_unfiltered = len(gene.polyA_internal_binding_sites)

    # read the bam
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    iterator = bam_file.fetch(gene.Chromosome, gene.Start, gene.End)
    int_prim_offsets = [] #store the position of the bam file for internal priming
    normal_prim_offsets = [] #store the position of the bam file for normal priming
    skipped_reads_offsets = [] #store the position of the bam file for skipped reads
    while True:
        try:
            offset = bam_file.tell() # store the current position of the bam file
            read = next(iterator)
            int_prim_flag = False
            if read.is_unmapped or read.is_supplementary:
                skipped_reads_offsets.append(offset)
                continue

            for start, end in gene.filtered_polyA_internal_binding_sites:
                if gene.Strand == '-':
                    # read.reference_start is 0-based inclusive — matches our 0-based site intervals
                    if start <= read.reference_start <= end:
                        int_prim_flag = True
                        break
                    elif read.reference_start < start:
                        break

                elif gene.Strand=='+':
                    if start <= read.reference_end <= end:
                        int_prim_flag = True
                        break
                    elif end < read.reference_end:
                        break
        

            if int_prim_flag:
                int_prim_offsets.append(offset) 
            else:
                normal_prim_offsets.append(offset)

    
        except StopIteration:
            break
    
    bam_file.close()

    # GET read counts as the pd.series
    rst_series = pd.Series({
        "Gene_id": gene.Gene_id,
        "IP_count": len(int_prim_offsets),
        "nonIP_count": len(normal_prim_offsets),
        "Skipped_reads": len(skipped_reads_offsets),
        "n_site_filtered": n_site_filtered,
        "n_site_unfiltered": n_site_unfiltered
    })

    return int_prim_offsets, normal_prim_offsets, skipped_reads_offsets, rst_series

def yield_read_from_offsets(bam_file_path: str, offsets: int):
    """
    Yield reads from a BAM file at the given byte offsets.
    Seeks to each offset individually and reads one record — this avoids the
    fragile pattern of mixing a shared iterator with seek(), which can return
    the wrong record depending on pysam's internal state.
    """
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    f = bam_file.fetch()
    for idx, i in enumerate(sorted(offsets)):
        bam_file.seek(i)
        try:
            read = next(f)
        except StopIteration:
            print(f"Warning: read not found at offset {i}")
            continue
        yield read

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
        with open(output_file, 'a') as f:
            f.write(sam_buffer)

def main():
    args = parsers.parse_arg()
    bam_template = pysam.AlignmentFile(args.bam_file, "rb")
    bam_header = bam_template.text
    # Warning message and exit if the output file already exists
    if args.output_sam:
        if helper.check_file(args.output_sam):
            helper.warning_msg(f"Output file {args.output_sam} already exists. Please remove it or use a different name.", printit=True)
            sys.exit(1)

    write_buffer(bam_header, args.output_sam)

    # Load optional gene subset
    gene_subset = None
    if args.gene_list:
        with open(args.gene_list) as f:
            gene_subset = set(line.strip() for line in f if line.strip())
        helper.green_msg(f"Loaded {len(gene_subset)} Gene_id(s) from {args.gene_list}.", printit=True)

    # Expect counts in the counter:
        # - internal_priming
        # - normal_priming
        # - skipped_reads
    # rst_counter = Counter()

    # Start
    gene_iter = GtfGeneClassGenerator(args.gtf_file, gene_subset=gene_subset)
    int_prim_offsets = []
    normal_prim_offsets = []
    skipped_reads_offsets = []
    gene_level_summary = []
    
     # Single process
    if args.processes == 1:
        for gene in tqdm(gene_iter):
            
            int_prim_offsets_gene, normal_prim_offsets_gene, skipped_reads_offsets_gene, rst_series = \
                process_gene(gene, args.bam_file, args.genome_ref, flank_buffer=args.flank_buffer)
            int_prim_offsets.extend(int_prim_offsets_gene)
            normal_prim_offsets.extend(normal_prim_offsets_gene)
            skipped_reads_offsets.extend(skipped_reads_offsets_gene)
            gene_level_summary.append(rst_series)
        # remove duplicates 
        int_prim_offsets = set(int_prim_offsets)
        normal_prim_offsets = set(normal_prim_offsets) - int_prim_offsets
        skipped_reads_offsets = set(skipped_reads_offsets) - int_prim_offsets - normal_prim_offsets
        
        rst_counter = {
            "internal_priming": len(int_prim_offsets),
            "normal_priming": len(normal_prim_offsets),
            "skipped_reads": len(skipped_reads_offsets)
        }
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
                            fasta_file_path=args.genome_ref,
                            flank_buffer=args.flank_buffer):
        
            int_prim_offsets_gene, normal_prim_offsets_gene, skipped_reads_offsets_gene,  rst_series= \
                future.result()
            
            int_prim_offsets.extend(int_prim_offsets_gene)
            normal_prim_offsets.extend(normal_prim_offsets_gene)
            skipped_reads_offsets.extend(skipped_reads_offsets_gene)
            gene_level_summary.append(rst_series)
            
        # remove duplicates 
        int_prim_offsets = set(int_prim_offsets)
        normal_prim_offsets = set(normal_prim_offsets) - int_prim_offsets
        skipped_reads_offsets = set(skipped_reads_offsets) - int_prim_offsets - normal_prim_offsets
        
        rst_counter = {
            "internal_priming": len(int_prim_offsets),
            "normal_priming": len(normal_prim_offsets),
            "skipped_reads": len(skipped_reads_offsets)
        }

    
    print(rst_counter)
    # Write the output sam file
    sam_tag = args.sam_tag
    # write the internal priming reads
    for read in yield_read_from_offsets(args.bam_file, int_prim_offsets):
        read.set_tag(sam_tag, 'T', value_type='A')
        write_buffer(read.to_string() + '\n', args.output_sam)
    # write the normal priming reads
    for read in yield_read_from_offsets(args.bam_file, normal_prim_offsets):
        read.set_tag(sam_tag, 'F', value_type='A')
        write_buffer(read.to_string() + '\n', args.output_sam)
    # write the skipped reads
    for read in yield_read_from_offsets(args.bam_file, skipped_reads_offsets):
        read.set_tag(sam_tag, 'N', value_type='A')
        write_buffer(read.to_string() + '\n', args.output_sam)

    # output the global summary
    total_mapped = rst_counter['internal_priming'] + rst_counter['normal_priming']
    proportion = rst_counter['internal_priming'] / total_mapped if total_mapped > 0 else 0.0
    with open(args.output_summary, 'w') as f:
        f.write(
            f"Total number of reads (mapped): {total_mapped}\n"
            f"Internally primed reads: {rst_counter['internal_priming']}\n"
            f"Normally primed reads: {rst_counter['normal_priming']}\n"
            f"Proportion of internal priming reads: {proportion}\n"
        )
        
    
    # output the gene level count
    gene_level_summary_df = pd.DataFrame(gene_level_summary)
    gene_level_summary_df.to_csv(args.output_gene_count, sep="\t", index=False)


    # Finished
    print("Finished")
if __name__ == "__main__":
    # redirect the print function to the stderror
    stdout = sys.stdout
    sys.stdout = sys.stderr
    main()