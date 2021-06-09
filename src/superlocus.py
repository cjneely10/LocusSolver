from typing import Dict, List


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


class SuperLocus:
    def __init__(self):
        self._start = -1
        self._end = -1
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
        if feature.start < self._start:
            self._start = feature.start
        if feature.end > self._end:
            self._end = feature.end
