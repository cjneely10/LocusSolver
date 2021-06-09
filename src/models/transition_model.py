from typing import List

from src.models.annotation_model import AnnotationModel
from src.util.superlocus_list import SuperLocusList


class TransitionModel:
    def __init__(self, annotation_models: List[AnnotationModel]):
        self._models = annotation_models
        self._merge_to_superloci()

    def _merge_to_superloci(self):
        self.feature_slls = {}
        for model in self._models:
            for (record_id, feature_list) in model.to_features():
                if record_id not in self.feature_slls.keys():
                    self.feature_slls[record_id] = feature_list
                else:
                    self.feature_slls[record_id].extend(feature_list)
        for record_id, feature_list in self.feature_slls.items():
            self.feature_slls[record_id] = SuperLocusList(feature_list)
