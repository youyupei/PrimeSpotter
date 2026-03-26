import os
import argparse
import textwrap

import helper

# arg parser
def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(helper.bold_text(
        '''PrimeSpotter: analyse and investigate internal priming in long-read sequencing.\n
        ''')),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--bam_file', type=ParserChecker.input_fn, help='Path to the BAM file', required=True)
    parser.add_argument('--gtf_file', type=ParserChecker.input_fn, help='Path to the GTF/GFF file', required=True)
    parser.add_argument('--genome-ref', type=ParserChecker.input_fn, help='Path to the FASTA file for genome reference',required=True)
    parser.add_argument('--processes', type=ParserChecker.processes, default=1, help='Number of processes to use')
    # add a output file prefix
    parser.add_argument('--output-sam', type=str, default="", help='Output filename, the output sam file will be write to the stdout if not specified')
    parser.add_argument('--output-summary', type=str, default="internal_priming_summary.txt", help='Output summary file')
    parser.add_argument('--output-gene-count', type=str, default="internal_priming_gene_count.tsv", help='Count of internal priming sites and reads per gene')
    parser.add_argument('--sam-tag', type=ParserChecker.sam_tag, default="IP", help='SAM tag for the whether it is an internal priming site read')
    parser.add_argument('--flank-buffer', type=int, default=10,
                        help='Number of bases to expand each detected poly-A site on both sides. '
                             'Larger values make the classifier more permissive (more reads flagged as IP). '
                             'Reduce to 2-5 if IP site locations look too broad in IGV.')
    parser.add_argument('--gene-list', type=ParserChecker.input_fn, default=None,
                        help='Optional path to a plain text file with one gene per line. '
                             'Accepts gene names (e.g. TTN, KCNQ1OT1) or Ensembl gene IDs '
                             '(e.g. ENSG00000155657), or a mix of both. '
                             'ENST transcript IDs are NOT supported — use gene names or ENSG IDs. '
                             'When provided, only these genes will be analysed.')

    args = parser.parse_args()
    return args

class ParserChecker:
    def sam_tag(tag):
        """
        Check if a given string is a valid SAM tag.

        Parameters:
        tag (str): The tag to check.

        Returns:
        str: The tag if valid.

        Raises:
        argparse.ArgumentTypeError: If the tag is not valid.
        """
        if isinstance(tag, str) and len(tag) == 2 and tag.isalpha() and tag.isupper():
            return tag
        else:
            raise argparse.ArgumentTypeError(f"Invalid SAM tag: '{tag}'. A valid SAM tag must be exactly 2 uppercase letters.")
    def processes(n):
        n = int(n)
        if n <= 0:
            raise argparse.ArgumentTypeError(f"Number of processes should be a positive integer, but got {n}")
        return n
    def input_fn(fn):
        if not os.path.isfile(fn):
            raise argparse.ArgumentTypeError(f"File '{fn}' does not exist.")
        return fn