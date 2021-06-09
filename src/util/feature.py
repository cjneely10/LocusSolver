class Feature:
    def __init__(self, start: int, end: int, strand: int):
        self._start = start
        self._end = end
        self._strand = strand

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    def __str__(self):
        return f"<Feature start: {self.start} end: {self.end} strand: {self.strand}>"

    def __repr__(self):
        return str(self)