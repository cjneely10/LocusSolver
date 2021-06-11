from typing import List, Optional


class Feature:
    def __init__(self, start: int, end: int, strand: int, _id: str, exon_list: Optional[List] = None):
        if exon_list is None:
            self._exons = []
        else:
            self._exons = exon_list
        self._start = start
        self._end = end
        self._strand = strand
        self._id = _id

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def id(self):
        return self._id

    @property
    def exons(self):
        return self._exons

    def __str__(self):
        return f"<Feature start: {self.start} end: {self.end} strand: {self.strand} sub-features: {self.exons}>"

    def __repr__(self):
        return str(self)
