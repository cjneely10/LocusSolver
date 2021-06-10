from typing import List, Dict, Callable, Generator

from Bio.SeqRecord import SeqRecord

from src.models.annotation_model import AnnotationModel
from src.util.superlocus import SuperLocus


class LocusFilter:
    def __init__(self, annotation_models: List[AnnotationModel]):
        self._models = annotation_models
        self.feature_slls: Dict[str, List[SuperLocus]] = AnnotationModel.merge(annotation_models)
        # for contig_id, sll in self.feature_slls.items():
        #     print(contig_id)
        #     for super_locus in sll:
        #         print(super_locus)

    def filter(self, fxn: Callable[[SuperLocus], SeqRecord]) -> Generator[SeqRecord, None, None]:
        pass
