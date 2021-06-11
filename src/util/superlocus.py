from sys import maxsize
from typing import Dict, List, Optional

from Bio.SeqFeature import SeqFeature


class SuperLocus:
    def __init__(self, feature: Optional[SeqFeature] = None, identifier: Optional[str] = None):
        self._start = maxsize
        self._end = -1 * maxsize
        # Store features by strand
        self._features: Dict[int, Dict[str, List[SeqFeature]]] = {1: {}, -1: {}}
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

    def add_feature(self, feature: SeqFeature, identifier: str):
        strand_dict = self._features[feature.strand]
        if identifier not in strand_dict.keys():
            strand_dict[identifier] = []
        strand_dict[identifier].append(feature)
        if feature.location.start < self._start:
            self._start = feature.location.start
        if feature.location.end > self._end:
            self._end = feature.location.end

    def overlaps(self, feature: SeqFeature) -> bool:
        return min(self.end, feature.location.end) > max(self.start, feature.location.start)

    def __gt__(self, other: SeqFeature) -> bool:
        return self.start > other.location.end

    def __lt__(self, other: SeqFeature) -> bool:
        return self.end < other.location.start

    def __str__(self):
        return f"<SuperLocus start: {self.start} end: {self.end} " \
               f"n_features(-): " \
               f"{','.join([name + '=' +  str(len(features)) for name, features in self.features[-1].items()])}; " \
               f"n_features(+): " \
               f"{','.join([name + '=' +  str(len(features)) for name, features in self.features[1].items()])}>"

    def __repr__(self):
        return str(self)
