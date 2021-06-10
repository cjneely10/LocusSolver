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

    def filter(self, fxn: Callable[[str, SuperLocus], FilterResult]) -> Generator[SeqRecord, None, None]:
        for i, (contig_id, sll) in enumerate(self._feature_slls.items(), start=1):
            record: SeqRecord = Annotation.genome_dict[contig_id]
            record.features = []
            for j, super_locus in enumerate(sll, start=1):
                result = fxn(contig_id, super_locus)
                for (identifier, selected_features) in result:
                    for k, selected_feature in enumerate(selected_features, start=1):
                        qualifiers = {
                            "source": identifier,
                            "ID": f"gene{i}.{j}.{k}-{identifier}-{contig_id}"
                                  f"-start:{selected_feature.start + 1}-end:{selected_feature.end}",
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

            yield record
