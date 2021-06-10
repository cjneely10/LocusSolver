from typing import Dict

import numpy as np

from src.util.annotation import Annotation

from src.filters import FilterResult
from src.util.feature import Feature
from src.util.superlocus import SuperLocus


def longest_cds(contig_id: str, super_locus: SuperLocus) -> FilterResult:
    feature: Feature
    out: FilterResult = []
    for strand in (1, -1):
        max_cds = (0, ("", []))  # length, FilterResult
        for name, features in super_locus.features[strand].items():
            total_cds_length = 0
            for feature in features:
                for exon in feature.exons:
                    total_cds_length += (exon.location.end - exon.location.start + 1)
            if total_cds_length > max_cds[0]:
                max_cds = (total_cds_length, (name, features))
        if len(max_cds[1][1]) > 0:
            out.append(max_cds[1])
    return out


def longest_range(contig_id: str, super_locus: SuperLocus) -> FilterResult:
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


def highest_scoring_intron_model(contig_id: str, super_locus: SuperLocus) -> FilterResult:
    feature: Feature
    out: FilterResult = []
    for strand in (1, -1):
        best_model = (-100000000000, ("", []))  # model_probability, FilterResult
        for name, features in super_locus.features[strand].items():
            total_score = 0.0
            for feature in features:
                score = 0.0
                for i in range(len(feature.exons) - 1):
                    exon1 = feature.exons[i]
                    exon2 = feature.exons[i + 1]
                    q_string = str(Annotation.subset_seq(
                            Annotation.genome_dict[contig_id], (exon1.location.end, exon2.location.start))).upper()
                    # if len(q_string) == Annotation.front_size + Annotation.end_size:
                    score += np.log10(Annotation.stored_models[name].intron_model.score(q_string))
                total_score += score
            if 10**total_score > best_model[0]:
                best_model = (10**total_score, (name, features))
        if len(best_model[1][1]) > 0:
            out.append(best_model[1])
    return out


class PriorityFilter:
    def __init__(self, priority_mapping: Dict[str, int]):
        self._priority_mapping = priority_mapping

    def __call__(self, contig_id: str, super_locus: SuperLocus) -> FilterResult:
        pass
