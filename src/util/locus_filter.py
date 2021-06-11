from typing import List, Dict, Callable, Generator

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from src.filters import FilterResult
from src.util.annotation import Annotation
from src.util.superlocus import SuperLocus


class LocusFilter:
    def __init__(self, annotation_models: List[Annotation]):
        self._models = annotation_models
        self._feature_slls: Dict[str, List[SuperLocus]] = Annotation.merge(annotation_models)

    def filter(self, fxn: Callable[[SuperLocus], FilterResult]) -> Generator[SeqRecord, None, None]:
        i = 1
        for (contig_id, sll) in self._feature_slls.items():
            record: SeqRecord = Annotation.genome_dict[contig_id]
            record.features = []
            for super_locus in sll:
                result = fxn(super_locus)
                for (identifier, selected_features) in result:
                    for selected_feature in selected_features:
                        qualifiers = {
                            "source": identifier,
                            "ID": f"gene{i}-{identifier}-{str(selected_feature.id)}"
                        }
                        sub_qualifiers = {"source": identifier}
                        top_feature = SeqFeature(
                            FeatureLocation(selected_feature.start, selected_feature.end),
                            type="gene",
                            strand=selected_feature.strand,
                            qualifiers=qualifiers
                        )
                        top_feature.sub_features = []
                        for sub_feature in selected_feature.exons:
                            top_feature.sub_features.append(
                                SeqFeature(
                                    FeatureLocation(sub_feature.location.start, sub_feature.location.end),
                                    type=sub_feature.type,
                                    strand=sub_feature.strand,
                                    qualifiers=sub_qualifiers)
                            )
                        record.features.append(top_feature)
                        i += 1

            yield record
