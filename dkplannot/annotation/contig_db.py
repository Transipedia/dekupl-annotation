import json
import sys
import os

class ContigDB(object):

    def __init__(self, tmp_dir: str) -> None:
        self.tmp_dir = tmp_dir
        # Create dir if it does not exist yet
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

    def getContig(self, tag: str):
        "Get a contig with its tag, if not contigs available create one"
        tag_file = self._getFilename(tag)
        contig_info = []
        if os.path.isfile(tag_file):
            f = open(tag_file, "r")
            json_line = f.readline()
            contig_info = json.loads(json_line)
        return contig_info

    def saveContig(self, tag: str, contig_info):
        tag_file = self._getFilename(tag)
        f = open(tag_file, "w+")
        f.write(json.dumps(contig_info))
    
    def _getFilename(self, tag: str):
        return self.tmp_dir + "/" + tag
