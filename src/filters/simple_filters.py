from typing import Dict

from Bio.SeqFeature import SeqFeature

from src.filters import FilterResult
from src.util.superlocus import SuperLocus


def longest_cds(super_locus: SuperLocus) -> FilterResult:
    feature: SeqFeature
    out: FilterResult = []
    for strand in (1, -1):
        max_cds = (0, [])  # length, FilterResult
        for name, features in super_locus.features[strand].items():
            total_cds_length = 0
            for feature in features:
                for exon in feature.sub_features:
                    total_cds_length += (exon.location.end - exon.location.start + 1)
            if total_cds_length > max_cds[0]:
                max_cds = (total_cds_length, features)
        if len(max_cds[1][1]) > 0:
            out.append(max_cds[1])
    return out


def longest_range(super_locus: SuperLocus) -> FilterResult:
    feature: SeqFeature
    out: FilterResult = []
    for strand in (1, -1):
        max_range = (0, [])  # length, FilterResult
        for name, features in super_locus.features[strand].items():
            longest = features[-1].end - features[0].start + 1
            if longest > max_range[0]:
                max_range = (longest, features)
        if len(max_range[1][1]) > 0:
            out.append(max_range[1])
    return out


class PriorityFilter:
    def __init__(self, priority_mapping: Dict[int, str]):
        self._priority_mapping = priority_mapping

    def __call__(self, super_locus: SuperLocus) -> FilterResult:
        feature: SeqFeature
        out: FilterResult = []
        for strand in (1, -1):
            priority = {name: features for name, features in super_locus.features[strand].items()}
            for i in range(1, len(self._priority_mapping) + 1):
                features = priority.get(self._priority_mapping[i])
                if features is not None:
                    out.append(features)
                    break
        return out


class Tier(PriorityFilter):
    def __init__(self, tier: int, priority_mapping: Dict[int, str]):
        super().__init__(priority_mapping)
        self._tier = tier

    def __call__(self, super_locus: SuperLocus) -> FilterResult:
        feature: SeqFeature
        out: FilterResult = []
        for strand in (1, -1):
            priority = {name: features for name, features in super_locus.features[strand].items()}
            if len(priority) < self._tier:
                continue
            for i in range(1, len(self._priority_mapping) + 1):
                features = priority.get(self._priority_mapping[i])
                if features is not None:
                    out.append(features)
                    break
        return out
