import argparse
import sys
import os

#from annotation import merged_counts
from annotation.merged_counts import MergedCounts
from annotation.gsnap import Gsnap

# Override the parser to print the help message in case of missing args
class DkplParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# Return the parser with all argment defined
def get_parser():
    parser = DkplParser(description='Annotation of dekupl-contigs')
    # Remove the original group and add 2 args group (required and optional)
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument(
        "-a, --contigs", 
        dest = "contigs",
        metavar = 'FILE',
        type = str,
        required = True,
        help = "merged-diff-counts.tsv.gz (contigs in \"{A}_vs_{B}_kmer_counts\" directory from Dekupl-run result)"
    )
    required.add_argument(
        "-o, --output", 
        dest = "output_directory",
        metavar = 'DIR',
        type = str,
        required = True,
        help = "Path to the output directory"
    )
    required.add_argument(
        "-i, --index", 
        dest = "index_directory",
        metavar = 'DIR',
        type = str,
        required = True,
        help = "Path to the index directory (created with dkplannot index)"
    )
    optional.add_argument(
        "-h, --help",
        help = "show this help message and exit"
    )
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()

    index_dir = "toy/index"
    index_name = "GRCh38-chr22"

    gsnap_index_dir = index_dir + "/gsnap"

    merged_counts = MergedCounts(args.contigs)
    gsnap = Gsnap(gsnap_index_dir, index_name)

    # Create output directory
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    else:
        sys.stderr.write("Output directory already exists")

    # Generate fasta contigs
    contigs_fasta = args.output_directory + "/contigs.fa.gz"
    if not os.path.isfile(contigs_fasta):
        merged_counts.generate_fasta(contigs_fasta)
    
    # Generate BAM with GSNAP
    contigs_bam = args.output_directory + "/contigs.bam"
    gsnap.generate_bam(contigs_fasta,contigs_bam)
    


main()