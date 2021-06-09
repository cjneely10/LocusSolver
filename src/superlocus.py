from typing import Dict, List


class Feature:
    def __init__(self, start: int, end: int, strand: int):
        self.start = start
        self.end = end
        self.strand = strand


class SuperLocus:
    def __init__(self, start: int, end: int):
        self._start = start
        self._end = end
        # Store features by strand
        self._features: Dict[int, Dict[str, List[Feature]]] = {1: {}, -1: {}}

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
