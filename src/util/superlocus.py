from sys import maxsize
from typing import Dict, List, Optional

from src.util.feature import Feature


class SuperLocus:
    def __init__(self, feature: Optional[Feature] = None, identifier: Optional[str] = None):
        self._start = maxsize
        self._end = -1 * maxsize
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

    @property
    def features(self):
        return self._features

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

    def overlaps(self, feature: Feature) -> bool:
        return min(self.end, feature.end) > max(self.start, feature.start)

    def __gt__(self, other: Feature) -> bool:
        return self.start > other.end

    def __lt__(self, other: Feature) -> bool:
        return self.end < other.start

    def __str__(self):
        return f"<SuperLocus start: {self.start} end: {self.end} " \
               f"n_features(-): " \
               f"{','.join([name + '=' +  str(len(features)) for name, features in self.features[-1].items()])}; " \
               f"n_features(+): " \
               f"{','.join([name + '=' +  str(len(features)) for name, features in self.features[1].items()])}>"

    def __repr__(self):
        return str(self)
