from typing import List

from src.models.annotation_model import AnnotationModel
from src.util.superlocus_list import SuperLocusList


class TransitionModel:
    def __init__(self, annotation_models: List[AnnotationModel]):
        self._models = annotation_models
        self._merge_to_superloci()

    def _merge_to_superloci(self):
        for model in self._models:
            for (record_id, feature_list) in model.to_features():
                sll = SuperLocusList(feature_list)
                print(record_id)
                for super_locus in sll.sorted:
                    print(super_locus)
