from typing import List, Dict, Callable, Generator, Tuple

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from src.util.annotation import Annotation
from src.util.feature import Feature
from src.util.superlocus import SuperLocus


class LocusFilter:
    def __init__(self, annotation_models: List[Annotation]):
        self._models = annotation_models
        self.feature_slls: Dict[str, List[SuperLocus]] = Annotation.merge(annotation_models)

    def filter(self, fxn: Callable[[SuperLocus], Tuple[str, Feature]]) -> Generator[SeqRecord, None, None]:
        for contig_id, sll in self.feature_slls.items():
            record: SeqRecord = Annotation.genome_dict[contig_id]
            record.features = []
            for i, super_locus in enumerate(sll):
                identifier, selected_feature = fxn(super_locus)
                qualifiers = {
                    "source": identifier,
                    "ID": f"gene{i}",
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
                            FeatureLocation(sub_feature.start, sub_feature.end),
                            type=sub_feature.type,
                            strand=sub_feature.strand,
                            qualifiers=sub_qualifiers)
                    )
                record.features.append(top_feature)

            yield record
