import sys
from typing import Dict, List, Optional


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


class SuperLocus:
    def __init__(self, feature: Optional[Feature] = None, identifier: Optional[str] = None):
        self._start = sys.maxsize
        self._end = -1 * sys.maxsize
        # Store features by strand
        self._features: Dict[int, Dict[str, List[Feature]]] = {1: {}, -1: {}}
        if feature is not None and identifier is not None:
            self.add_feature(feature, identifier)

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    def add_feature(self, feature: Feature, identifier: str):
        strand_dict = self._features[feature.strand]
        if identifier not in strand_dict.keys():
            strand_dict[identifier] = []
        strand_dict[identifier].append(feature)
        strand_dict[identifier].sort(key=lambda feat: feat.start)
        if feature.start < self._start:
            self._start = feature.start
        if feature.end > self._end:
            self._end = feature.end

    def overlaps(self, feature: Feature):
        return min(self.end, feature.end) > max(self.start, feature.start)

    def __gt__(self, other: Feature):
        return self.start > other.end

    def __lt__(self, other: Feature):
        return self.end < other.start

    def __str__(self):
        return f"<SuperLocus start: {self.start} end: {self.end} features: {str(self._features)}>"

    def __repr__(self):
        return str(self)
