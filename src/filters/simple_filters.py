from typing import Dict

from src.filters import FilterResult
from src.util.feature import Feature
from src.util.superlocus import SuperLocus


def longest_cds(super_locus: SuperLocus) -> FilterResult:
    feature: Feature
    out: FilterResult = []
    for strand in (1, -1):
        max_cds = (0, ("", []))  # length, FilterResult
        for name, features in super_locus.features[strand].items():
            total_cds_length = 0
            for feature in features:
                total_cds_length += (feature.end - feature.start + 1)
            if total_cds_length > max_cds[0]:
                max_cds = (total_cds_length, (name, features))
        if len(max_cds[1][1]) > 0:
            out.append(max_cds[1])
    return out


def longest_range(super_locus: SuperLocus) -> FilterResult:
    feature: Feature
    out: FilterResult = []
    for strand in (1, -1):
        max_range = (0, ("", []))  # length, FilterResult
        for name, features in super_locus.features[strand].items():
            longest = features[-1].end - features[0].start + 1
            if longest > max_range[0]:
                max_range = (longest, (name, features))
        if len(max_range[1][1]) > 0:
            out.append(max_range[1])
    return out


def highest_scoring_intron_model(super_locus: SuperLocus) -> FilterResult:
    pass


class PriorityFilter:
    def __init__(self, priority_mapping: Dict[str, int]):
        self._priority_mapping = priority_mapping

    def __call__(self, super_locus: SuperLocus) -> FilterResult:
        pass
