import os

from annotation.merged_counts import MergedCounts

MATRIX_TOY = "tests/datas/merged-diff-counts.tsv.gz"
TEMP_DIR = "tests/datas/tmp"

def test_mergedcounts_class():
    output_file = TEMP_DIR + "/merged_counts.fa.gz"
    merged_counts = MergedCounts(MATRIX_TOY)
    if not os.path.isdir(TEMP_DIR):
        os.mkdir(TEMP_DIR)

    elif os.path.exists(output_file):
        os.remove(output_file)

    merged_counts.generate_fasta(output_file)
    assert os.path.isfile(output_file)
