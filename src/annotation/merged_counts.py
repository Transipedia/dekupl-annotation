import gzip

class MergedCounts(object):

    def __init__(self, file: str) -> None:
        self._file = file

    @property
    def file(self) -> str:
        return self._file

    @file.setter
    def file(self, value: str) -> None:
        if value is not None:
            self._file = value

        else:
            raise ValueError("Can't set a NoneType value.")

    def generate_fasta(self, output: str) -> None:
        """
        Generate a fasta from dekupl output matrix
        """
        out = gzip.open(output, 'wt')
        with gzip.open(self.file, "rt") as content:
            header = content.readline()
            for line in content:
                split_line = line.rstrip().split("\t")
                out.write(">" + split_line[2] + "\n")
                out.write(split_line[1] + "\n")

        out.close()
