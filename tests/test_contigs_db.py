import os

from annotation.contig_db import ContigDB

tmp_dir = ".contig_db"

contig_db = ContigDB(tmp_dir)

contig = {
    'toto' : 'tata'
}

contig2 = {
    'toto' : 'titi'
}

# Write contig
contig_db.saveContig("AAA",contig)

# Overload contig
contig_db.saveContig("AAA",contig2)

# Load contig
contig_loaded = contig_db.getContig("AAA")

# Verify data
assert(contig_loaded['toto'] == "titi")