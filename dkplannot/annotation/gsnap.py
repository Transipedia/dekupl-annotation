import sys
import subprocess

class Gsnap(object):

    def __init__(self, index_dir: str, index_name: str, threads = 1) -> None:
        self.index_dir = index_dir
        self.index_name = index_name
        self.threads = threads
    
    def generate_bam(self, fasta : str, output: str) -> None:
        """
        Generate alignments in BAM format from a FASTA input file
        """
        command = " ".join([
            "gsnap -t ", str(self.threads),
            "-A sam", # -A, --format=STRING            Another format type, other than default.
            "-N 1", # -N, --novelsplicing=INT              Look for novel splicing (0=no (default), 1=yes)
            "-D", self.index_dir, #  -D, --dir=directory            Genome directory.
            "-d", self.index_name, # -d, --db=STRING                Genome database
            #"-w 50000", # -w, --localsplicedist=INT            Definition of local novel splicing event (default 200000)
            "--gunzip", # Uncompress gzipped input files
            fasta,
            "| samtools view -bh >", # Convert output SAM to BAM format with samtools
            output
        ])

        sys.stderr.write("Executing GSNAP\nCommand: %s\n" % command)
        subprocess.call(command,shell=True)
        